#! /usr/bin/env python3

import json, os, tempfile, unittest

from bam_qc_metrics import bam_qc, fast_metric_finder

class test(unittest.TestCase):

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
        self.reference = None
        self.n_as_mismatch = False
        self.maxDiff = None # uncomment to show the (very long) full output diff

    def test_default_analysis(self):
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path, self.markdup_path,
                    self.n_as_mismatch, self.quality, self.reference, sample_rate=None, tmpdir=self.tmpdir, verbose=False)
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
            self.assertEqual(expected_variables[key], output[key])
        # now check all output data
        with (open(self.expected_path)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()
        
    def test_downsampled_analysis(self):
        sample_rate = 10
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path,
                    self.markdup_path, self.n_as_mismatch, self.quality, self.reference,
                    sample_rate, tmpdir=self.tmpdir, verbose=False)
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
            self.assertEqual(expected_variables[key], output[key])
        # now check all output data
        with (open(self.expected_path_downsampled)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()

    def test_missing_inputs(self):
        # test possible missing inputs:
        # - ESTIMATED_LIBRARY_SIZE in mark duplicates text
        # - FFQ/LFQ in samtools stats
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path,
                    self.markdup_path, self.n_as_mismatch, self.quality, self.reference,
                    sample_rate=None, tmpdir=self.tmpdir, verbose=False)
        # for low-coverage runs, ESTIMATED_LIBRARY_SIZE value is missing from mark duplicates text
        # test input file also has variant '## METRICS CLASS ...' line
        metrics_found = qc.read_mark_duplicates_metrics(self.markdup_path_low_cover)
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
        fast_finder = fast_metric_finder(self.bam_path, self.reference, self.insert_max)
        fq_result = fast_finder.fq_stats([])
        fq_expected = ({},{})
        self.assertEqual(fq_expected, fq_result)
        qc.cleanup()
        
    def tearDown(self):
        self.tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
