#!/usr/bin/env python

"""
# =============================================================================

Tests for the k-mer cache (neptune.KMerCache) and its use by Neptune.

# =============================================================================
"""

import os
import sys
import shutil
import unittest

from tests.TestingUtility import *

import neptune.KMerCache as KMerCache
import neptune.CountKMers as CountKMers
import neptune.ExtractSignatures as ExtractSignatures
import neptune.JobManagerParallel as JobManagerParallel
import neptune.Neptune as Neptune
import neptune.Execution as Execution

OUTPUT = getPath("tests/output/kmer_cache/temp.dir")
SIMPLE = getPath("tests/data/neptune/simple.fasta")
ALTERNATIVE = getPath("tests/data/neptune/alternative.fasta")


def buildJobManager():

    os.makedirs(OUTPUT, exist_ok=True)
    return JobManagerParallel.JobManagerParallel(OUTPUT, OUTPUT, 2)


def readAll(baseLocation, organization):

    contents = ""
    for location in KMerCache.fileLocations(baseLocation, organization):
        with open(location, 'r') as f:
            contents += f.read()
    return contents


class CountingJobManager():
    """Wraps a job manager, recording which genomes were counted."""

    def __init__(self, jobManager):
        self.jobManager = jobManager
        self.counted = []

    def createCountJob(self, inputLocation, outputLocation, k, organization):
        self.counted.append(inputLocation)
        return self.jobManager.createCountJob(
            inputLocation, outputLocation, k, organization)

    def runJobs(self, jobs):
        self.jobManager.runJobs(jobs)


class TestPopulate(unittest.TestCase):

    def tearDown(self):
        shutil.rmtree(OUTPUT, ignore_errors=True)

    def test_populate_matches_count(self):

        for organization in [0, 1]:

            cache = os.path.join(OUTPUT, "cache")
            locations = KMerCache.populate(
                buildJobManager(), cache, 5, organization, [SIMPLE])

            expected = os.path.join(OUTPUT, "expected.kmers")
            CountKMers.count(SIMPLE, expected, 5, organization)

            self.assertEqual(
                readAll(locations[0], organization),
                readAll(expected, organization))
            self.assertTrue(
                KMerCache.isCached(cache, 5, organization, SIMPLE))

            # no temporary files left behind
            directory = KMerCache.cacheDirectory(cache, 5, organization)
            self.assertFalse(any(
                name.startswith(KMerCache.TEMP_PREFIX)
                for name in os.listdir(directory)))

    def test_reuse(self):

        cache = os.path.join(OUTPUT, "cache")
        manager = CountingJobManager(buildJobManager())

        KMerCache.populate(manager, cache, 5, 1, [SIMPLE])
        self.assertEqual(manager.counted, [SIMPLE])

        manager.counted = []
        KMerCache.populate(manager, cache, 5, 1, [SIMPLE, ALTERNATIVE])
        self.assertEqual(manager.counted, [ALTERNATIVE])

        # a different k is a different cache entry
        manager.counted = []
        KMerCache.populate(manager, cache, 6, 1, [SIMPLE])
        self.assertEqual(manager.counted, [SIMPLE])

    def test_stale_entry_recounted(self):

        cache = os.path.join(OUTPUT, "cache")
        genome = os.path.join(OUTPUT, "genome.fasta")
        os.makedirs(OUTPUT, exist_ok=True)
        shutil.copy(SIMPLE, genome)

        manager = CountingJobManager(buildJobManager())
        KMerCache.populate(manager, cache, 5, 0, [genome])

        with open(genome, 'a') as f:
            f.write("GGGGGGGGGG\n")

        self.assertFalse(KMerCache.isCached(cache, 5, 0, genome))

        manager.counted = []
        locations = KMerCache.populate(manager, cache, 5, 0, [genome])
        self.assertEqual(manager.counted, [genome])

        expected = os.path.join(OUTPUT, "expected.kmers")
        CountKMers.count(genome, expected, 5, 0)
        self.assertEqual(readAll(locations[0], 0), readAll(expected, 0))

    def test_duplicate_names(self):

        other = os.path.join(OUTPUT, "other", "simple.fasta")
        os.makedirs(os.path.dirname(other))
        shutil.copy(ALTERNATIVE, other)

        with self.assertRaises(RuntimeError):
            KMerCache.populate(
                buildJobManager(), os.path.join(OUTPUT, "cache"), 5, 0,
                [SIMPLE, other])


class TestNeptuneWithCache(unittest.TestCase):
    """Counting + aggregation via the cache matches the uncached result."""

    def tearDown(self):
        shutil.rmtree(OUTPUT, ignore_errors=True)

    def aggregate(self, output, cache, organization):

        parameters = {
            ExtractSignatures.INCLUSION: [SIMPLE],
            ExtractSignatures.EXCLUSION: [ALTERNATIVE],
            CountKMers.KMER: 5,
            CountKMers.ORGANIZATION: organization,
            KMerCache.CACHE: cache,
            Neptune.OUTPUT: output}

        execution = Execution.Execution(buildJobManager(), parameters)
        inclusion, exclusion = Neptune.countKMers(execution)
        Neptune.aggregateKMers(execution, inclusion, exclusion)

        with open(execution.aggregateLocation, 'r') as f:
            return f.read()

    def test_aggregate_identical(self):

        for organization in [0, 1]:

            cache = os.path.join(OUTPUT, "cache")

            uncached = self.aggregate(
                os.path.join(OUTPUT, "a"), None, organization)
            cached = self.aggregate(
                os.path.join(OUTPUT, "b"), cache, organization)
            reused = self.aggregate(
                os.path.join(OUTPUT, "c"), cache, organization)

            self.assertTrue(len(uncached) > 0)
            self.assertEqual(uncached, cached)
            self.assertEqual(uncached, reused)

            # cache survives aggregation
            self.assertTrue(
                KMerCache.isCached(cache, 5, organization, SIMPLE))
            self.assertTrue(
                KMerCache.isCached(cache, 5, organization, ALTERNATIVE))


class TestCacheCommand(unittest.TestCase):

    def tearDown(self):
        shutil.rmtree(OUTPUT, ignore_errors=True)

    def test_cache_subcommand(self):

        cache = os.path.join(OUTPUT, "cache")
        sys.argv[1:] = [
            "cache", "-i", SIMPLE, ALTERNATIVE, "-k", "5",
            "--organization", "1", "--kmer-cache", cache, "-p", "2"]

        Neptune.main()

        self.assertTrue(KMerCache.isCached(cache, 5, 1, SIMPLE))
        self.assertTrue(KMerCache.isCached(cache, 5, 1, ALTERNATIVE))

    def test_cache_subcommand_requires_k(self):

        sys.argv[1:] = [
            "cache", "-i", SIMPLE, "--kmer-cache",
            os.path.join(OUTPUT, "cache")]

        with self.assertRaises(SystemExit):
            Neptune.main()


if __name__ == '__main__':
    unittest.main()
