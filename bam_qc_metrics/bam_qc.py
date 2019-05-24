#! /usr/bin/env python3

"""Main script and class to compute BAM QC metrics"""

import argparse, csv, json, re, pysam, sys

class bam_qc:

    METADATA_KEYS = [
        "barcode",
        "instrument",
        "lane",
        "library",
        "run name",
        "sample"
    ]
    
    def __init__(self, bam_path, metadata_path, ref_path):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        self.bam_path = bam_path
        self.ref_path = ref_path
        self.samtools_metrics = self.run_samtools()
        print(self.samtools_metrics)

    def run_samtools(self):
        result = pysam.stats(self.bam_path)
        # parse output from 'samtools stats', filtering out comment lines
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", result)),
            delimiter="\t"
        )
        # read some simple metrics
        keys = ['RL', 'FRL', 'LRL']
        metrics = {}
        for row in reader:
            if row[0] in keys:
                metrics[row[0]] = [int(n) for n in row[1:]]
        return metrics
        
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
    parser.add_argument('-b', '--bam', metavar='PATH',
                        help='Path to input BAM file')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata')
    parser.add_argument('-o', '--out', metavar='PATH',
                        help='Path for JSON output, or - for STDOUT')
    parser.add_argument('-r', '--reference', metavar='PATH',
                        help='Path to reference BED file')
    args = parser.parse_args()
    qc = bam_qc(args.bam, args.metadata, args.reference)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
