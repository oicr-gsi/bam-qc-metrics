#! /usr/bin/env python3

import json, os, tempfile, unittest

from bam_qc_metrics import bam_qc

class test(unittest.TestCase):

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory(prefix='bam_qc_test')
        self.tmpdir = self.tmp.name
        self.testdir = os.path.join(os.path.dirname(__file__))
        self.metadata_path = os.path.join(self.testdir, 'metadata.json')
        self.bam_path = os.path.join(self.testdir, 'neat_5x_EX_hg19_chr21.bam')
        self.markdup_path = os.path.join(self.testdir, 'marked_dup_metrics.txt')
        self.target_path = os.path.join(self.testdir,'SureSelect_All_Exon_V4_Covered_Sorted_chr21.bed')
        expected_path = os.path.join(self.testdir, 'expected.json')
        with (open(expected_path)) as f: self.expected = json.loads(f.read())
        
    def test(self):
        qc = bam_qc(self.bam_path, self.target_path, self.metadata_path, self.markdup_path)
        out_path = os.path.join(self.tmpdir, 'out.json')
        qc.write_output(out_path)
        self.assertTrue(os.path.exists(out_path))
        with (open(out_path)) as f: output = json.loads(f.read())
        self.assertEqual(output, self.expected)

    def tearDown(self):
        self.tmp.cleanup()


if __name__ == '__main__':
    unittest.main()
