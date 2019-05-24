#! /usr/bin/env python3

"""Main script and class to compute BAM QC metrics"""

import argparse, json, sys

class bam_qc:

    METADATA_KEYS = [
        "barcode",
        "instrument",
        "lane",
        "library",
        "run name",
        "sample"
    ]
    
    def __init__(self, metadata_path):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        
    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        out_file = None
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata')
    parser.add_argument('-o', '--out', metavar='PATH',
                        help='Path for JSON output, or - for STDOUT')
    args = parser.parse_args()
    qc = bam_qc(args.metadata)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
