"""
   Work with windowed coverage from mosdepth, produce a coverage histogram
   should be used only with call-ready (multilane) metrics as the original code
   handles lane-level metrics well
"""
import pandas as pd
import json
import re
import argparse


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Run parsing script to generate assay scan report table')
    parser.add_argument('-f', '--files', help='Comma separated files with coverage data', required=True)
    parser.add_argument('-o', '--out', help='Output', required=False, default="test_coverage_histogram.json")
    args = parser.parse_args()

    input_files = re.split(',', args.files)
    output_file = args.out
    """
       Read all data, append to array of data pandas frames
    """
    dfs = []
    for f in input_files:
        dfs.append(pd.read_csv(f, sep="\t", header=None))

    """
       Massage the data back into 4-column DataFrame summing data of each 4th column
       and using this sum as a value in the col #4 of the resulting DataFrame
    """
    concat_dfs = pd.concat(dfs, axis=1)
    summed_df = dfs[0].iloc[:, :3].copy()
    summed_df.columns = ['chrom', 'start', 'stop']
    value_cols = [df.iloc[:, 3] for df in dfs]
    summed_df['coverage'] = pd.concat(value_cols, axis=1).sum(axis=1)
    summed_df.drop(summed_df[summed_df.coverage == 0].index, inplace=True)
    summed_df.reset_index(inplace=True)
    """
        Our goal is to calculate how many bases we have at a certain depth of coverage
        we calculate rounded coverage (int) so we need to quickly create a column
        we need max rounded coverage with non-zero bases covered
    """
    summed_df['rounded_coverage'] = round(summed_df['coverage'])
    histogram = {}
    window = summed_df['stop'][0] - summed_df['start'][0]
    """
       Count the rows covered at certain depth to speed things up
       This is faster than for loop! We will use it for safe retrieval
       Count the rows covered at certain depth to speed things up, max_coverage is the largest coverage 
       with at least 3 windows covered at this coverage (otherwise the data are too large and unlikely 
       to be very useful towards the right tail)
       
    """
    base_counts = summed_df["rounded_coverage"].value_counts().to_dict()
    filtered_coverages = filter(lambda item: item[1] >= 3, base_counts.items())
    f_coverages = dict(filtered_coverages)
    max_coverage = round(max(list(f_coverages.keys())))

    for depth in range(0, max_coverage + 1):
        rows_covered = base_counts.get(depth, 0)
        histogram[depth] = int(window * rows_covered)
    """
        Write it all down into a json file. This can directly go into bamQC as the
        format is identical to what the original histogram-building code produces    
    """
    out = open(output_file, "w")
    json.dump(histogram, out, sort_keys=True)
    out.close()
