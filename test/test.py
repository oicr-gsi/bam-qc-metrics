#! /usr/bin/env python3

import json, os, tempfile, unittest

from bam_qc_metrics import bam_qc

class test(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='bam_qc_test_')
        self.tmpdir = self.tmp.name
        self.testdir = os.path.join(os.path.dirname(__file__))
        self.metadata_path = os.path.join(self.testdir, 'metadata.json')
        self.bam_path = os.path.join(self.testdir, 'neat_5x_EX_hg19_chr21.bam')
        self.markdup_path = os.path.join(self.testdir, 'marked_dup_metrics.txt')
        self.target_path = os.path.join(self.testdir,'SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed')
        self.expected_path = os.path.join(self.testdir, 'expected.json')
        self.expected_path_downsampled = os.path.join(self.testdir, 'expected_downsampled.json')
        #self.maxDiff = None # uncomment to show the (very long) full output diff
        
    def test(self):
        quality = 30
        qc = bam_qc(self.bam_path, self.target_path, self.metadata_path, self.markdup_path, quality)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        with (open(self.expected_path)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.tmp.cleanup()
        
    def test_downsample(self):
        quality = 30
        sample_rate = 10
        qc = bam_qc(self.bam_path, self.target_path, self.metadata_path, self.markdup_path, quality,
                    sample_rate=sample_rate)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        with (open(self.expected_path_downsampled)) as f: expected = json.loads(f.read())
        self.assertEqual(output, expected)
        qc.tmp.cleanup()
        
    def tearDown(self):
        self.tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
