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
    
    def __init__(self, bam_path, target_path, metadata_path=None, mark_duplicates_path=None):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        self.bam_path = bam_path
        self.target_path = target_path
        self.samtools_metrics = self.run_samtools()
        self.mark_duplicates_metrics = {
            "ESTIMATED_LIBRARY_SIZE": None,
            "HISTOGRAM": {},
            "LIBRARY": None,
            "PERCENT_DUPLICATION": None,
            "READ_PAIRS_EXAMINED": None,
            "READ_PAIR_DUPLICATES": None,
            "READ_PAIR_OPTICAL_DUPLICATES": None,
            "UNMAPPED_READS": None,
            "UNPAIRED_READS_EXAMINED": None,
            "UNPAIRED_READ_DUPLICATES": None
        }
        if mark_duplicates_path != None:
            self.mark_duplicates_metrics = self.read_mark_duplicates_metrics(mark_duplicates_path)

    def read_mark_duplicates_metrics(self, input_path):
        section = 0
        line_count = 0
        with open(input_path) as f: lines = f.readlines()
        keys = []
        values = []
        hist = {}
        for line in lines:
            line_count += 1
            line = line.strip()
            if re.match('## METRICS CLASS\s+net\.sf\.picard\.sam\.DuplicationMetrics', line):
                section += 1
            elif section == 1:
                keys = re.split("\t", line)
                section += 1
            elif section == 2:
                values = re.split("\t", line)
                section += 1
            elif section == 3 and re.match('## HISTOGRAM', line):
                section += 1
            elif section == 4 and re.match("BIN\tVALUE", line):
                section += 1
            elif section == 5 and re.match("[0-9]+\.[0-9]+\t[0-9]+\.{0,1}[0-9]*$", line):
                terms = re.split("\t", line)
                # JSON doesn't allow numeric dictionary keys, so hist_bin is stringified in output
                hist_bin = round(float(terms[0])) if re.search('\.0$', terms[0]) else terms[0]
                hist[hist_bin] = float(terms[1])
            elif re.match('#', line) and section < 4 or line == '':
                continue
            else:
                params = (input_path, section, line)
                msg = "Failed to parse duplicate metrics path %s, section %d, line %d" % params
                raise ValueError(msg)
        if len(keys) != len(values):
            raise ValueError("Key and value lists from %s are of unequal length" % input_path)
        metrics = {}
        for i in range(len(keys)):
            if keys[i] == 'PERCENT_DUPLICATION':
                metrics[keys[i]] = float(values[i])
            elif keys[i] == 'LIBRARY':
                metrics[keys[i]] = values[i]
            else:
                metrics[keys[i]] = int(values[i])
        metrics['HISTOGRAM'] = hist
        return metrics

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
            'reads unmapped': 'unmapped reads',
            'non-primary alignments': 'non primary reads',
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
        output['mark duplicates'] = self.mark_duplicates_metrics
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Picard MarkDuplicates. Optional.')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-o', '--out', metavar='PATH',
                        help='Path for JSON output, or - for STDOUT')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional; if given, must be sorted in same order as BAM file.')
    args = parser.parse_args()
    qc = bam_qc(args.bam, args.target, args.metadata, args.mark_duplicates)
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
