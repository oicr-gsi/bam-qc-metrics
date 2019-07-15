#! /usr/bin/env python3

import json, os, tempfile, unittest

from bam_qc_metrics import bam_qc

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
        #self.maxDiff = None # uncomment to show the (very long) full output diff

    def test(self):
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path,
                    self.markdup_path, self.quality)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        with (open(self.expected_path)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()
        
    def test_downsample(self):
        sample_rate = 10
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path,
                    self.markdup_path, self.quality, sample_rate=sample_rate)
        out_path = os.path.join(self.tmpdir, 'out_downsampled.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        with (open(self.expected_path_downsampled)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.cleanup()

    def test_missing_library_size(self):
        # for low-coverage runs, ESTIMATED_LIBRARY_SIZE value is missing from mark duplicates text
        # test input file also has variant '## METRICS CLASS ...' line
        qc = bam_qc(self.bam_path, self.target_path, self.insert_max, self.metadata_path,
                    self.markdup_path_low_cover, self.quality)
        metrics_found = qc.read_mark_duplicates_metrics(self.markdup_path_low_cover)
        with (open(self.expected_metrics_low_cover)) as f: metrics_expected = json.loads(f.read())
        # Don't directly compare found/expected HISTOGRAM values. Keys are integers and strings, respectively.
        # (Annoyingly, JSON format insists dictionary keys must be strings)
        histogram_found = metrics_found['HISTOGRAM']
        histogram_expected = metrics_expected['HISTOGRAM']
        self.assertEqual(len(histogram_found), len(histogram_expected))
        for key in histogram_found.keys():
            self.assertEqual(histogram_found[key], histogram_expected[str(key)])
        del metrics_found['HISTOGRAM']
        del metrics_expected['HISTOGRAM']
        self.assertEqual(metrics_found, metrics_expected)
        qc.cleanup()
        
    def tearDown(self):
        self.tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
