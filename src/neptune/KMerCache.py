#!/usr/bin/env python

"""
# =============================================================================

Copyright Government of Canada 2015-2024

Written by: Eric Marinier, Public Health Agency of Canada,
    National Microbiology Laboratory

Funded by the National Micriobiology Laboratory and the Genome Canada / Alberta
    Innovates Bio Solutions project "Listeria Detection and Surveillance
    using Next Generation Genomics"

Licensed under the Apache License, Version 2.0 (the "License"); you may not use
this file except in compliance with the License. You may obtain a copy of the
License at:

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed
under the License is distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR
CONDITIONS OF ANY KIND, either express or implied. See the License for the
specific language governing permissions and limitations under the License.

# =============================================================================
"""

"""
# =============================================================================

K-MER CACHE

This module manages a persistent cache of per-genome k-mer count files, so
that k-mers counted for a genome in one Neptune execution can be reused in
later executions, regardless of whether the genome is used as an inclusion or
an exclusion target.

A genome's k-mer files depend only on the genome, the k-mer size and the
degree of organization. Cache entries are therefore stored as:

[CACHE]/k[K]_o[ORGANIZATION]/[GENOME BASENAME].kmers[.TAG]
[CACHE]/k[K]_o[ORGANIZATION]/[GENOME BASENAME].kmers.done

The ".done" marker is written only after all of a genome's k-mer files are
complete. It records the absolute path, size and modification time of the
source genome; an entry whose marker does not match the current genome file
is considered stale and is recounted.

Entries are first written under a temporary name and renamed into place, so
concurrent Neptune executions sharing a cache never observe partial files.

# =============================================================================
"""

import argparse
import json
import os
import uuid

import neptune.Utility as Utility
import neptune.CountKMers as CountKMers

"""
# =============================================================================

GLOBALS

# =============================================================================
"""

KMERS_EXTENSION = ".kmers"
DONE_EXTENSION = ".done"
TEMP_PREFIX = ".tmp-"

PROGRAM_DESCRIPTION = "Counts the k-mers of one or more genomes and stores \
    them in a k-mer cache directory. The cache can then be passed to Neptune \
    with --kmer-cache, so that the genomes' k-mers are not recounted in every \
    Neptune execution. Genomes already present in the cache are skipped."

PROGRAM_USAGE = "%(prog)s cache -i INPUT [INPUT ...] -k K --kmer-cache CACHE"

LONG = "--"
SHORT = "-"

# Cache directory
CACHE = "kmer_cache"
CACHE_LONG = LONG + "kmer-cache"
CACHE_HELP = "The directory of the k-mer cache. Per-genome k-mer counts are \
    read from and written to this directory, so that they may be reused \
    between Neptune executions. The k-mer size (-k) should be specified when \
    using a cache; otherwise it is estimated from the inclusion genomes and \
    may differ between executions, preventing reuse."

# Input genomes (cache command)
INPUT = "input"
INPUT_LONG = LONG + INPUT
INPUT_SHORT = SHORT + "i"
INPUT_HELP = "The FASTA files, or directories of FASTA files, to count \
    and add to the cache."

# Parallelization (cache command)
PARALLELIZATION = "parallelization"
PARALLELIZATION_LONG = LONG + PARALLELIZATION
PARALLELIZATION_SHORT = SHORT + "p"
PARALLELIZATION_HELP = "The number of processes to run simultaneously."

ORGANIZATION_DEFAULT = 3
PARALLELIZATION_DEFAULT = 8

"""
# =============================================================================

CACHE DIRECTORY

Returns the directory holding cache entries for the given k and organization.

# =============================================================================
"""
def cacheDirectory(cacheLocation, k, organization):

    return os.path.abspath(
        os.path.join(
            cacheLocation, "k" + str(k) + "_o" + str(organization)))


"""
# =============================================================================

CACHE LOCATION

Returns the k-mer file base location of the given genome within the cache.
This location is what is passed to the aggregation step, which appends the
organization tags itself when organization > 0.

# =============================================================================
"""
def cacheLocation(cacheRoot, k, organization, genomeLocation):

    return os.path.join(
        cacheDirectory(cacheRoot, k, organization),
        os.path.basename(genomeLocation) + KMERS_EXTENSION)


"""
# =============================================================================

FILE LOCATIONS

Returns every file that makes up a k-mer entry with the given base location.

# =============================================================================
"""
def fileLocations(baseLocation, organization):

    if organization == 0:
        return [baseLocation]

    return [
        baseLocation + "." + tag
        for tag in Utility.getAggregationTags(organization)]


"""
# =============================================================================

SOURCE SIGNATURE

Describes the genome file, so that stale cache entries can be detected.

# =============================================================================
"""
def sourceSignature(genomeLocation):

    status = os.stat(genomeLocation)

    return {
        "source": os.path.abspath(genomeLocation),
        "size": status.st_size,
        "mtime_ns": status.st_mtime_ns}


"""
# =============================================================================

IS CACHED

Returns whether a complete and current cache entry exists for the genome.

# =============================================================================
"""
def isCached(cacheRoot, k, organization, genomeLocation):

    baseLocation = cacheLocation(cacheRoot, k, organization, genomeLocation)
    doneLocation = baseLocation + DONE_EXTENSION

    if not os.path.isfile(doneLocation):
        return False

    try:
        with open(doneLocation, 'r') as doneFile:
            recorded = json.load(doneFile)

    except (OSError, ValueError):
        return False

    if recorded != sourceSignature(genomeLocation):
        return False

    return all(
        os.path.isfile(location)
        for location in fileLocations(baseLocation, organization))


"""
# =============================================================================

CHECK UNIQUE NAMES

Cache entries are keyed by the genome's file name, so two different genomes
with the same file name cannot share a cache.

# =============================================================================
"""
def checkUniqueNames(genomeLocations):

    seen = {}

    for location in genomeLocations:

        name = os.path.basename(location)
        location = os.path.abspath(location)

        if name in seen and seen[name] != location:
            raise RuntimeError(
                "ERROR: The k-mer cache requires genome file names to be "
                "unique, but '" + name + "' refers to both " + seen[name]
                + " and " + location + ".\n")

        seen[name] = location


"""
# =============================================================================

POPULATE

Counts the k-mers of every genome that is not already cached and adds them to
the cache. Returns the cache base locations of all requested genomes, in the
order they were given.

[JOB MANAGER] [jobManager]
    Used to run the CountKMers jobs in parallel.

# =============================================================================
"""
def populate(jobManager, cacheRoot, k, organization, genomeLocations):

    checkUniqueNames(genomeLocations)

    directory = cacheDirectory(cacheRoot, k, organization)
    os.makedirs(directory, exist_ok=True)

    jobs = []
    pending = []    # (genome, temporary base, final base)
    locations = []
    scheduled = set()

    for genomeLocation in genomeLocations:

        finalLocation = cacheLocation(
            cacheRoot, k, organization, genomeLocation)
        locations.append(finalLocation)

        if finalLocation in scheduled or isCached(
                cacheRoot, k, organization, genomeLocation):
            continue

        scheduled.add(finalLocation)

        temporaryLocation = os.path.join(
            directory,
            TEMP_PREFIX + uuid.uuid4().hex + "-"
            + os.path.basename(finalLocation))

        jobs.append(jobManager.createCountJob(
            genomeLocation, temporaryLocation, k, organization))
        pending.append((genomeLocation, temporaryLocation, finalLocation))

    print("k-mer cache: " + str(len(genomeLocations) - len(pending))
          + " cached, " + str(len(pending)) + " to count.")

    try:
        jobManager.runJobs(jobs)

    except BaseException:
        for _, temporaryLocation, _ in pending:
            removeEntry(temporaryLocation, organization)
        raise

    for genomeLocation, temporaryLocation, finalLocation in pending:
        commit(genomeLocation, temporaryLocation, finalLocation, organization)

    return locations


"""
# =============================================================================

COMMIT

Moves a finished temporary entry into place and writes its ".done" marker
last, so that the entry only becomes visible once it is complete.

# =============================================================================
"""
def commit(genomeLocation, temporaryLocation, finalLocation, organization):

    # Invalidate any previous (stale) entry before replacing its files.
    doneLocation = finalLocation + DONE_EXTENSION
    if os.path.exists(doneLocation):
        os.remove(doneLocation)

    for temporary, final in zip(
            fileLocations(temporaryLocation, organization),
            fileLocations(finalLocation, organization)):
        os.replace(temporary, final)

    temporaryDone = temporaryLocation + DONE_EXTENSION
    with open(temporaryDone, 'w') as doneFile:
        json.dump(sourceSignature(genomeLocation), doneFile)

    os.replace(temporaryDone, doneLocation)


"""
# =============================================================================

REMOVE ENTRY

Removes any files of a (temporary) entry.

# =============================================================================
"""
def removeEntry(baseLocation, organization):

    for location in fileLocations(baseLocation, organization) \
            + [baseLocation + DONE_EXTENSION]:
        if os.path.exists(location):
            os.remove(location)


"""
# =============================================================================

MAIN (neptune cache ...)

# =============================================================================
"""
def main(arguments, version):

    import neptune.JobManagerParallel as JobManagerParallel

    parser = argparse.ArgumentParser(
        prog="neptune cache",
        description=PROGRAM_DESCRIPTION)

    parser.add_argument(
        SHORT + "V", LONG + "version",
        action='version',
        version='%(prog)s ' + str(version))

    required = parser.add_argument_group("REQUIRED")

    required.add_argument(
        INPUT_SHORT, INPUT_LONG,
        dest=INPUT, help=INPUT_HELP,
        type=str, required=True, nargs='+')

    required.add_argument(
        CountKMers.KMER_SHORT, CountKMers.KMER_LONG,
        dest=CountKMers.KMER, help=CountKMers.KMER_HELP,
        type=int, required=True)

    required.add_argument(
        CACHE_LONG,
        dest=CACHE, help=CACHE_HELP,
        type=str, required=True)

    optional = parser.add_argument_group("OPTIONAL")

    optional.add_argument(
        CountKMers.ORGANIZATION_LONG,
        dest=CountKMers.ORGANIZATION,
        help=CountKMers.ORGANIZATION_HELP
        + " Must match the value used by the Neptune executions that use "
        "this cache.",
        type=int, default=ORGANIZATION_DEFAULT)

    optional.add_argument(
        PARALLELIZATION_SHORT, PARALLELIZATION_LONG,
        dest=PARALLELIZATION, help=PARALLELIZATION_HELP,
        type=int, default=PARALLELIZATION_DEFAULT)

    args = parser.parse_args(arguments)

    if args.kmer < 1:
        parser.error("The k-mer size is out of range.")

    if args.organization < 0:
        parser.error("The organization is out of range.")

    genomeLocations = []
    Utility.expandInput(args.input, genomeLocations)

    cacheRoot = os.path.abspath(args.kmer_cache)
    logLocation = os.path.join(cacheRoot, "log")
    os.makedirs(logLocation, exist_ok=True)

    jobManager = JobManagerParallel.JobManagerParallel(
        cacheRoot, logLocation, args.parallelization)

    populate(
        jobManager, cacheRoot, args.kmer, args.organization, genomeLocations)

    print("k-mer cache ready: "
          + cacheDirectory(cacheRoot, args.kmer, args.organization))
