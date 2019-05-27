#! /usr/bin/env python3

"""Main script and class to compute BAM QC metrics"""

import argparse, csv, json, re, pysam, sys

class bam_qc:

    METADATA_KEYS = [
        'barcode',
        'instrument',
        'lane',
        'library',
        'run name',
        'sample'
    ]
    
    def __init__(self, bam_path, metadata_path, ref_path):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        self.bam_path = bam_path
        self.ref_path = ref_path
        self.samtools_metrics = self.run_samtools()

    def run_samtools(self):
        result = pysam.stats(self.bam_path)
        reader = csv.reader(
            filter(lambda line: line!="" and line[0]!='#', re.split("\n", result)),
            delimiter="\t"
        )
        # summary numbers (SN) fields denoted in float_keys are floats; integers otherwise
        float_keys = [
            'error rate',
            'average quality',
            'insert size average',
            'insert size standard deviation',
            'percentage of properly paired reads (%)'
        ]
        # map from SN field names to output keys
        key_map = {
            'bases mapped (cigar)': 'bases mapped',
            'average length': 'average read length',
            'insert size average': 'insert mean',
            'insert size standard deviation': 'insert stdev',
            'reads mapped': 'mapped reads',
            'reads mapped and paired': 'reads mapped and paired',
            'mismatches': 'mismatched bases',
            # 'reads paired' > 0 => 'paired end' = True
            'reads paired': 'paired reads',
            'pairs on different chromosomes': 'pairsMappedToDifferentChr',
            'reads properly paired': 'properly paired reads',
            'raw total sequences': 'total reads',
            'reads unmapped': 'unmapped reads'
        }
        samtools_stats = {}
        samtools_stats['inserted bases'] = 0
        samtools_stats['deleted bases'] = 0
        read_len = {'total': 0, 'count': 0}
        for row in reader:
            if row[0] == 'SN':
                samtools_key = re.sub(':$', '', row[1])
                if samtools_key not in key_map: continue
                if samtools_key in float_keys: val = float(row[2])
                else: val = int(row[2])
                samtools_stats[key_map[samtools_key]] = val
            elif row[0] == 'ID':
                samtools_stats['inserted bases'] += int(row[1]) * int(row[2])
                samtools_stats['deleted bases'] += int(row[1]) * int(row[3])
            elif row[0] == 'RL':
                length = int(row[1])
                count = int(row[2])
                read_len['total'] += length * count
                read_len['count'] += count
        samtools_stats['paired end'] = samtools_stats['paired reads'] > 0
        if read_len['count'] > 0:
            samtools_stats['average read length'] = float(read_len['total']) / read_len['count']
        else:
            samtools_stats['average read length'] = None
        return samtools_stats
        
    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        for key in self.samtools_metrics.keys():
            output[key] = self.samtools_metrics.get(key)
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
