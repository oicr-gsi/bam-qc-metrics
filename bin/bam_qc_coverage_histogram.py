"""
   Work with summary and global distribution data mosdepth, produce a coverage histogram
   should be used only with lane-level (single lane) metrics (the original code)
"""
import json
import argparse
import csv

"""
   Original code from the bamQC WDL
"""
def make_histogram(summary_file_path: str, global_dist_file_path: str, out_file: str):
    summary = open(summary_file_path, 'r', newline='')
    globalDist = open(global_dist_file_path, 'r', newline='')

    # read chromosome lengths from the summary
    summaryReader = csv.reader(summary, delimiter="\t")
    lengthByChr = {}

    for row in summaryReader:
        if row[0] == 'chrom' or row[0] == 'total':
            continue # skip initial header row, and final total row
        if row[0].endswith('_region'):
            continue # skip contigs from target file, if passed
        lengthByChr[row[0]] = int(row[1])
    chromosomes = sorted(lengthByChr.keys())
    summary.close()

    # read the cumulative distribution for each chromosome
    globalReader = csv.reader(globalDist, delimiter="\t")
    cumDist = {}
    for k in chromosomes:
        cumDist[k] = {}
    for row in globalReader:
        if row[0] == "total":
            continue
        cumDist[row[0]][int(row[1])] = float(row[2])
    globalDist.close()

    # convert the cumulative distributions to non-cumulative and populate histogram
    # if the input BAM is empty, chromosomes and histogram will also be empty
    histogram = {}
    for k in chromosomes:
        depths = sorted(cumDist[k].keys())
        dist = {}
        for i in range(len(depths)-1):
            depth = depths[i]
            nextDepth = depths[i+1]
            dist[depth] = cumDist[k][depth] - cumDist[k][nextDepth]
        maxDepth = max(depths)
        dist[maxDepth] = cumDist[k][maxDepth]
        # now find the number of loci at each depth of coverage to construct the histogram
        for depth in depths:
            loci = int(round(dist[depth]*lengthByChr[k], 0))
            histogram[depth] = histogram.get(depth, 0) + loci

    # if histogram is non-empty, fill in zero values for missing depths
    for i in range(max(histogram.keys(), default=0)):
        if i not in histogram:
            histogram[i] = 0
    out = open(out_file, "w")
    json.dump(histogram, out, sort_keys=True)
    out.close()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run this to construct coverage histogram for single-lane data')
    parser.add_argument('-s', '--summary', help='Coverage summary data', required=True)
    parser.add_argument('-g', '--gdist', help="Global distribution data", required=True)
    parser.add_argument('-o', '--out', help='Output', required=False, default="coverage_histogram.json")
    args = parser.parse_args()

    summary_file = args.summary
    global_dist_file = args.gdist
    output_file = args.out

    make_histogram(summary_file, global_dist_file, output_file)
