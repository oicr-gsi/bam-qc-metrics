"""
   This script is for handling bamQC report files (json format) for multiple lanes and merging them
   into a single report (call ready bamQC report)

   Some of the metrics need to be summed, some - averaged and some will go into the final report as-is.
   Only relevant metrics are kept, lane-specific metrics are not accepted. These metrics are lane-specific
   and cannot be used in call ready mode
"""

import json
import argparse
import sys
from re import split

DEFAULT_MERGED_REPORT = "bamQC_final_report.json"
ASIS = 0
MEAN = 1
SUM = 2
HIST = 3
MAX = 4

class bamqc_merger:
    PRECISION = 1
    """ Unknown and line-specific metrics will be skipped """
    supported_metrics = {
        "alignment reference": ASIS,
        "average read length": MEAN,
        "bases mapped": SUM,
        "bases per target": HIST,
        "coverage_histogram": HIST,
        "coverage per target": HIST,
        "deleted bases": SUM,
        "donor": ASIS,
        "downsampled total": SUM,
        "hard clip bases": SUM,
        "group id": ASIS,
        "insert max": MAX,
        "insert size average": MEAN,
        "insert size histogram": HIST,
        "insert size standard deviation": MEAN,
        "inserted bases": SUM,
        "instrument": ASIS,
        "library design": ASIS,
        "mapped reads": SUM,
        "mark duplicates": HIST,
        "mean genomic coverage": MEAN,
        "mean target coverage": MEAN,
        "mismatch by cycle": HIST,
        "mismatched bases": SUM,
        "non primary reads": SUM,
        "number of targets": MAX,
        "package version": ASIS,
        "paired end": ASIS,
        "paired reads": SUM,
        "pairsMappedAbnormallyFar": SUM,
        "pairsMappedToDifferentChr": SUM,
        "properly paired reads": SUM,
        "qual fail reads": SUM,
        "read 1 aligned by cycle": HIST,
        "read 1 average length": MEAN,
        "read 1 deletion by cycle": HIST,
        "read 1 hard clip by cycle": HIST,
        "read 1 insertion by cycle": HIST,
        "read 1 length histogram": HIST,
        "read 1 quality by cycle": HIST,
        "read 1 quality histogram": HIST,
        "read 1 soft clip by cycle": HIST,
        "read 2 aligned by cycle": HIST,
        "read 2 average length": MEAN,
        "read 2 deletion by cycle": HIST,
        "read 2 hard clip by cycle": HIST,
        "read 2 insertion by cycle": HIST,
        "read 2 length histogram": HIST,
        "read 2 quality by cycle": HIST,
        "read 2 quality histogram": HIST,
        "read 2 soft clip by cycle": HIST,
        "read ? aligned by cycle": HIST,
        "read ? average length": MEAN,
        "read ? deletion by cycle": HIST,
        "read ? hard clip by cycle": HIST,
        "read ? insertion by cycle": HIST,
        "read ? length histogram": HIST,
        "read ? quality by cycle": HIST,
        "read ? quality histogram": HIST,
        "read ? soft clip by cycle": HIST,
        "reads mapped and paired": SUM,
        "reads on target": SUM,
        "reads per start point": MEAN,
        "readsMissingMDtags": SUM,
        "soft clip bases": SUM,
        "target file": ASIS,
        "target sizes": HIST,
        "tissue origin": ASIS,
        "tissue type": ASIS,
        "total bases on target": SUM,
        "total reads": SUM,
        "total target size": MAX,
        "unmapped reads": SUM,
        "workflow version": ASIS
    }
    report = {}
    inputs = []

    ''' Function will accept a list of json files and return a list of dicts'''
    @staticmethod
    def load_inputs(report_files: list):
        replist = []
        for rf in report_files:
            with open(rf, "r") as repFile:
                replist.append(json.load(repFile))
            repFile.close()
        return replist

    ''' Data handling functions '''
    def _merge_as_is(self, metric: str, value, count: int, total_reads: int, running_total: int):
        self.report[metric] = value

    def _merge_mean(self, metric: str, value, count: int, total_reads: int, running_total: int):
        vetted_value = 0 if value is None else value
        if metric in self.report:
            current_normalized = self.report[metric] * running_total
            new_normalized = vetted_value * total_reads
            total_updated = running_total + total_reads
            '''
                Non-normalized implementation:
                self.report[metric] = round(self.report[metric] * (count-1)/count + vetted_value/count, PRECISION)
            '''
            self.report[metric] = round((current_normalized + new_normalized)/total_updated, self.PRECISION)
        else:
            self.report[metric] = vetted_value

    def _merge_sum(self, metric: str, value, count: int, total_reads: int, running_total: int):
        if value is not None and count > 0:
            self.report[metric] = self.report.get(metric, 0) + value

    def _merge_hist(self, metric: str, value, count: int, total_reads: int, running_total: int):
        """quality by cycle hist, cov. per target, targets sizes and PERCENT DUPLICATION are handled differently"""
        if not isinstance(value, dict):
            return
        if metric not in self.report:
            self.report[metric] = {}
        for entry, number in value.items():
            if metric.endswith('quality by cycle') or metric == "coverage per target" or entry == 'PERCENT_DUPLICATION':
                current_normalized = self.report[metric].get(entry, 0) * running_total
                new_normalized = number * total_reads
                total_updated = running_total + total_reads
                self.report[metric][entry] = round((current_normalized + new_normalized)/total_updated, self.PRECISION)
            elif metric == "target sizes":
                if self.report[metric].get(entry, 0) == 0:
                    self.report[metric][entry] = number
            else:
                self.report[metric][entry] = self.report[metric].get(entry, 0) + number

    def _merge_max(self, metric: str, value, count: int, total_reads: int, running_total: int):
        if metric in self.report.keys():
            self.report[metric] = max(self.report[metric], value)
        else:
            self.report[metric] = value

    '''Merging function'''
    def merge_reports(self):
        method_dispatch = {
            ASIS: self._merge_as_is,
            MEAN: self._merge_mean,
            SUM: self._merge_sum,
            HIST: self._merge_hist,
            MAX: self._merge_max
        }
        try:
            count = 0
            running_total_reads = 0
            for rep in self.inputs:
                count += 1
                total_reads = rep["total reads"]
                for metric, value in rep.items():
                    if metric not in self.supported_metrics:
                        continue
                    method = self.supported_metrics[metric]
                    handler = method_dispatch.get(method)
                    if handler:
                        handler(metric, value, count, total_reads, running_total_reads)
                running_total_reads += total_reads
        except Exception as e:
            sys.stderr.write("Error reading data: {0}".format(e))

    def get_report(self):
        return self.report

    def __init__(self, inputs: list):
        self.inputs = bamqc_merger.load_inputs(inputs)
        self.report = {}
        self.merge_reports()


def main():
    parser = argparse.ArgumentParser(description='QC metrics merger.')
    parser.add_argument('-l', '--list', help='List of JSON reports.', required=True)
    parser.add_argument('-o', '--output', help='', required=False, default=DEFAULT_MERGED_REPORT)
    args = parser.parse_args()

    report_list = args.list
    output = args.output

    reports = split(",", report_list)
    '''If we have only one input report, just dump it into the json'''
    if len(reports) == 1:
        report = reports[0]
    else:
        b_merger = bamqc_merger(reports)
        report = b_merger.get_report()

    with open(output, "w") as out:
        json.dump(report, out, ensure_ascii=True)


if __name__ == "__main__":
    main()
