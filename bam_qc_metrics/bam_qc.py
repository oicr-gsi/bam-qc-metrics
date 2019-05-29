#! /usr/bin/env python3

"""Main script and class to compute BAM QC metrics"""

import argparse, csv, json, os, re, pysam, sys

class bam_qc:

    METADATA_KEYS = [
        'barcode',
        'instrument',
        'lane',
        'library',
        'run name',
        'sample'
    ]
    
    def __init__(self, bam_path, target_path, metadata_path=None, mark_duplicates_path=None,
                 trim_quality=None):
        with open(metadata_path) as f: self.metadata = json.loads(f.read())
        self.bam_path = bam_path
        self.target_path = target_path
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
        self.trim_quality = trim_quality
        self.samtools_metrics = self.run_samtools()

    def fq_stats(self, rows, precision=6):
        '''
        Compute quality metrics from either FFQ or LFQ entries in samtools stats:
            - Mean quality by cycle
            - Quality histogram
        '''
        meanByCyc = {}
        max_width = max([len(row) for row in rows])
        histogram = {q: 0 for q in range(max_width-2)}
        for row in rows:
            cycle = int(row[1])
            counts = [int(n) for n in row[2:]]
            total = 0
            count = 0
            for qscore in range(len(counts)):
                total += counts[qscore]*qscore
                count += counts[qscore]
                histogram[qscore] += counts[qscore]
            meanByCyc[cycle] = round(float(total) / count, precision) if count > 0 else 0
        return (meanByCyc, histogram)

    def mean_read_length(self, rows):
        '''Process RL (read length), FRL (first RL) or LRL (last RL) rows for mean read length'''
        total = 0
        count = 0
        for row in rows:
            length = int(row[1])
            length_count = int(row[2])
            total += length * length_count
            count += length_count
        mean_rl = float(total) / count if count > 0 else 0
        return mean_rl

    def read_length_histogram(self, rows):
        '''Process RL (read length), FRL (first RL) or LRL (last RL) rows for read length histogram'''
        histogram = {}
        for row in rows:
            histogram[int(row[1])] = int(row[2])
        return histogram
    
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
                # but rounding removes the trailing '.0'
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
        if self.trim_quality == None:
            result = pysam.stats(self.bam_path)
        else:
            result = pysam.stats("-q", str(self.trim_quality), self.bam_path)
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
        samtools_stats['insert size histogram'] = {}
        read_len = {'total': 0, 'count': 0}
        ffq_rows = []
        lfq_rows = []
        rl_rows = []
        frl_rows = []
        lrl_rows = []
        for row in reader:
            if row[0] == 'FFQ':
                ffq_rows.append(row)
            elif row[0] == 'FRL':
                frl_rows.append(row)
            elif row[0] == 'ID':
                samtools_stats['inserted bases'] += int(row[1]) * int(row[2])
                samtools_stats['deleted bases'] += int(row[1]) * int(row[3])
            elif row[0] == 'IS':
                samtools_stats['insert size histogram'][int(row[1])] = int(row[2])
            elif row[0] == 'LFQ':
                lfq_rows.append(row)
            elif row[0] == 'LRL':
                lrl_rows.append(row)
            elif row[0] == 'RL':
                rl_rows.append(row)
            elif row[0] == 'SN':
                samtools_key = re.sub(':$', '', row[1])
                if samtools_key not in key_map: continue
                if samtools_key in float_keys: val = float(row[2])
                else: val = int(row[2])
                samtools_stats[key_map[samtools_key]] = val
        samtools_stats['average read length'] = self.mean_read_length(rl_rows)
        samtools_stats['paired end'] = samtools_stats['paired reads'] > 0
        samtools_stats['qual cut'] = self.trim_quality
        # TODO pysam.view may have to be adjusted for downsampling
        samtools_stats['qual fail reads'] = int(pysam.view('-c', self.bam_path)) \
                                            - samtools_stats['mapped reads']
        samtools_stats['read 1 average length'] = self.mean_read_length(frl_rows)
        samtools_stats['read 2 average length'] = self.mean_read_length(lrl_rows)
        samtools_stats['read 1 length histogram'] = self.read_length_histogram(frl_rows)
        samtools_stats['read 2 length histogram'] = self.read_length_histogram(lrl_rows)
        (ffq_mean_by_cycle, ffq_histogram) = self.fq_stats(ffq_rows)
        (lfq_mean_by_cycle, lfq_histogram) = self.fq_stats(lfq_rows)
        samtools_stats['read 1 quality by cycle'] = ffq_mean_by_cycle
        samtools_stats['read 1 quality histogram'] = ffq_histogram
        samtools_stats['read 2 quality by cycle'] = lfq_mean_by_cycle
        samtools_stats['read 2 quality histogram'] = lfq_histogram
        return samtools_stats
        
    def write_output(self, out_path):
        output = {}
        for key in self.METADATA_KEYS:
            output[key] = self.metadata.get(key)
        for key in self.samtools_metrics.keys():
            output[key] = self.samtools_metrics.get(key)
        output['mark duplicates'] = self.mark_duplicates_metrics
        output['target'] = os.path.abspath(self.target_path)
        output['target_size'] = None # TODO placeholder; find with pybedtools
        if out_path != '-':
            out_file = open(out_path, 'w')
        else:
            out_file = sys.stdout
        print(json.dumps(output), file=out_file)
        if out_path != '-':
            out_file.close()


def validate_args(args):
    valid = True
    q_threshold = None
    if args.trim_quality != None:
        try:
            q_threshold = int(args.trim_quality)
        except ValueError:
            sys.stderr.write("ERROR: Quality must be an integer.\n")
            valid = False
        if q_threshold < 0:
            sys.stderr.write("ERROR: Quality cannot be negative.\n")
            valid = False
    for path_arg in (args.bam, args.target, args.metadata, args.mark_duplicates):
        if path_arg == None:
            continue
        if not os.path.exists(path_arg):
            sys.stderr.write("ERROR: Path %s does not exist.\n" % path_arg)
            valid = False
        elif not os.path.isfile(path_arg):
            sys.stderr.write("ERROR: Path %s is not a file.\n" % path_arg)
            valid = False
        elif not os.access(path_arg, os.R_OK):
            sys.stderr.write("ERROR: Path %s is not readable.\n" % path_arg)
            valid = False
    if args.out != '-':
        # ugly but robust Python idiom to resolve path of parent directory
        parent_path = os.path.abspath(os.path.join(args.out, os.pardir))
        if not os.path.isdir(parent_path):
            sys.stderr.write("ERROR: Parent of %s is not a directory.\n" % args.out)
            valid = False
        elif not os.access(parent_path, os.W_OK):
            sys.stderr.write("ERROR: Parent directory of %s is not writable.\n" % args.out)
            valid = False
    return valid

def main():
    parser = argparse.ArgumentParser(description='QC for BAM files.')
    parser.add_argument('-b', '--bam', metavar='PATH', required=True,
                        help='Path to input BAM file. Required.')
    parser.add_argument('-d', '--mark-duplicates', metavar='PATH',
                        help='Path to text file output by Picard MarkDuplicates. Optional.')
    parser.add_argument('-m', '--metadata', metavar='PATH',
                        help='Path to JSON file containing metadata. Optional.')
    parser.add_argument('-o', '--out', metavar='PATH', required=True,
                        help='Path for JSON output, or - for STDOUT. Required.')
    parser.add_argument('-q', '--trim-quality', metavar='QSCORE',
                        help='Samtools threshold for trimming alignment quality. Optional.')
    parser.add_argument('-t', '--target', metavar='PATH',
                        help='Path to target BED file, containing targets to calculate coverage '+\
                        'against. Optional; if given, must be sorted in same order as BAM file.')
    args = parser.parse_args()
    if not validate_args(args): exit(1)
    qc = bam_qc(args.bam, args.target, args.metadata, args.mark_duplicates, int(args.trim_quality))
    qc.write_output(args.out)

if __name__ == "__main__":
    main()
