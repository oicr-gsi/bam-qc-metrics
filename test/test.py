#! /usr/bin/env python3

import json, os, sys, tempfile, unittest

from bam_qc_metrics import bam_qc, fast_metric_finder

class test(unittest.TestCase):

    OICR_REF = '/oicr/data/reference/genomes/homo_sapiens/UCSC/Genomic/UCSC_hg19/fasta/hg19.fa'
    LOCAL_REF = os.path.join(os.path.expanduser('~'), 'data/reference/hg19.fa')

    def setUp(self):
        self.quality = 30
        self.insert_max = 1500
        self.tmp = tempfile.TemporaryDirectory(prefix='bam_qc_test_')
        self.tmpdir = self.tmp.name
        self.testdir = os.path.dirname(os.path.realpath(__file__))
        self.datadir = os.path.realpath(os.path.join(self.testdir, '..', 'data'))
        self.metadata_path = os.path.join(self.datadir, 'metadata.json')
        self.bam_path = os.path.join(self.datadir, 'neat_5x_EX_hg19_chr21.bam')
        self.markdup_path = os.path.join(self.datadir, 'marked_dup_metrics.txt')
        self.markdup_path_low_cover = os.path.join(self.datadir, 'marked_dup_metrics_low_cover.txt')
        self.target_path = os.path.join(self.datadir,'SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed')
        self.expected_path = os.path.join(self.datadir, 'expected.json')
        self.expected_path_downsampled = os.path.join(self.datadir, 'expected_downsampled.json')
        self.expected_metrics_low_cover = os.path.join(self.datadir, 'expected_metrics_low_cover.json')
        if os.path.exists(self.OICR_REF):
            self.reference = self.OICR_REF
        elif os.path.exists(self.LOCAL_REF):
            self.reference = self.LOCAL_REF
        else:
            msg = "WARNING: Default reference file %s or %s not found" % \
                  (self.OICR_REF, self.LOCAL_REF)
            print(msg, file=sys.stderr)
            self.reference = None
        self.n_as_mismatch = False
        self.verbose = False
        self.maxDiff = None # uncomment to show the (very long) full output diff

    def test_default_analysis(self):
        config =  {
            "bam": self.bam_path,
            "target": self.target_path,
            "insert max": self.insert_max,
            "metadata": self.metadata_path,
            "mark duplicates": self.markdup_path,
            "n as mismatch": self.n_as_mismatch,
            "skip below mapq": self.quality,
            "reference": self.reference,
            "sample rate": None,
            "temp dir": self.tmpdir,
            "verbose": self.verbose
        }
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        # do individual sanity checks on some variables
        # helps validate results if expected output JSON file has been changed
        expected_variables = {
            "inserted bases": 315,
            "reads per start point": 1.031,
            "readsMissingMDtags": 80020,
            "sample rate": 1,
            "total reads": 80020,
            "total target size": 527189,
        }
        for key in expected_variables.keys():
            expected = expected_variables[key]
            got = output[key]
            try:
                self.assertEqual(expected, got)
            except AssertionError:
                print("\nFailed on metric '"+key+"': Expected", expected, ", got", got,
                      file=sys.stderr)
                raise
        # now check all output data
        with (open(self.expected_path)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()
        
    def test_downsampled_analysis(self):
        config =  {
            "bam": self.bam_path,
            "target": self.target_path,
            "insert max": self.insert_max,
            "metadata": self.metadata_path,
            "mark duplicates": self.markdup_path,
            "n as mismatch": self.n_as_mismatch,
            "skip below mapq": self.quality,
            "reference": self.reference,
            "sample rate": 10,
            "temp dir": self.tmpdir,
            "verbose": self.verbose
        }
        sample_rate = 10
        qc = bam_qc(config)
        out_path = os.path.join(self.tmpdir, 'out_downsampled.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        # do individual sanity checks on some variables
        # helps validate results if expected output JSON file has been changed
        expected_variables = {
            "inserted bases": 315,
            "reads per start point": 1.002, # downsampled
            "readsMissingMDtags": 7874, # downsampled
            "sample rate": 10,
            "total reads": 80020,
            "total target size": 527189,
        }
        for key in expected_variables.keys():
            expected = expected_variables[key]
            got = output[key]
            try:
                self.assertEqual(expected, got)
            except AssertionError:
                print("\nFailed on metric '"+key+"': Expected", expected, ", got", got,
                      file=sys.stderr)
                raise
        # now check all output data
        with (open(self.expected_path_downsampled)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()

    def test_missing_inputs(self):
        # test possible missing inputs:
        # - ESTIMATED_LIBRARY_SIZE in mark duplicates text
        # - FFQ/LFQ in samtools stats
        config =  {
            "bam": self.bam_path,
            "target": self.target_path,
            "insert max": self.insert_max,
            "metadata": self.metadata_path,
            "mark duplicates": self.markdup_path,
            "n as mismatch": self.n_as_mismatch,
            "skip below mapq": self.quality,
            "reference": self.reference,
            "sample rate": None,
            "temp dir": self.tmpdir,
            "verbose": self.verbose
        }
        qc = bam_qc(config)
        # for low-coverage runs, ESTIMATED_LIBRARY_SIZE value is missing from mark duplicates text
        # test input file also has variant '## METRICS CLASS ...' line
        metrics_found = qc.read_mark_dup(self.markdup_path_low_cover)
        with (open(self.expected_metrics_low_cover)) as f: metrics_expected = json.loads(f.read())
        # Found/expected HISTOGRAM keys are integers and strings, respectively.
        # (Annoyingly, JSON format insists dictionary keys must be strings)
        histogram_found = metrics_found['HISTOGRAM']
        histogram_expected = metrics_expected['HISTOGRAM']
        self.assertEqual(len(histogram_found), len(histogram_expected))
        for key in histogram_found.keys():
            self.assertEqual(histogram_found[key], histogram_expected[str(key)])
        del metrics_found['HISTOGRAM']
        del metrics_expected['HISTOGRAM']
        self.assertEqual(metrics_found, metrics_expected)
        # test empty FFQ/LFQ result from samtools stats; may occur for small input datasets
        # requires a fast_metric_finder object
        fast_finder = fast_metric_finder(self.bam_path,
                                         self.reference,
                                         self.insert_max,
                                         self.n_as_mismatch)
        fq_result = fast_finder.fq_stats([])
        fq_expected = ({},{})
        self.assertEqual(fq_expected, fq_result)
        qc.cleanup()
        
    def tearDown(self):
        self.tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
