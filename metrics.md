# Metric specifications for BamQC

BamQC metrics are computed using a number of methods, including third-party software tools, and output in JSON format. The output also contains metadata, such as the instrument and lane names; and input parameters, such as the target file. This document summarizes how metrics are computed, and details some special cases of interest.

Some metrics have been renamed in the move from the old Perl to new Python implementation.

## Summary of fields

| Name                            | Former name             | Source                  | DS |
| --------------------------------|-------------------------|-------------------------|----|
| average read length             |                         | samtools stats: RL      | N  |
| barcode                         |                         | metadata                | N  |
| bases mapped                    | aligned bases           | samtools stats: SN      | N  |
| deleted bases                   |                         | samtools stats: ID      | N  |
| hard clip bases                 |                         | CIGAR                   | Y  |
| insert max                      |                         | input parameter         | N  |
| insert size average             | insert mean             | samtools stats: SN      | N  |
| insert size histogram           | insert histogram        | samtools stats: IS      | N  |
| insert size standard deviation  | insert stdev            | samtools stats: SN      | N  |
| inserted bases                  |                         | samtools stats: ID      | N  |
| instrument                      |                         | metadata                | N  |
| lane                            |                         | metadata                | N  |
| library                         |                         | metadata                | N  |
| mapped reads                    |                         | samtools stats: SN      | N  |
| mark duplicates                 |                         | Picard MarkDuplicates   | N  |
| mismatched bases                | mismatch bases          | samtools stats: SN      | N  |
| non primary reads               |                         | samtools stats: SN      | N  |
| number of targets               |                         | bedtools                | N  |
| paired end                      | number of ends          | samtools stats: SN      | N  |
| paired reads                    |                         | samtools stats: SN      | N  |
| pairsMappedAbnormallyFar        |                         | insert size histogram   | Y  |
| pairsMappedToDifferentChr       |                         | samtools stats: SN      | N  |
| properly paired reads           |                         | samtools stats: SN      | N  |
| qual cut                        |                         | input parameter         | N  |
| qual fail reads                 |                         | samtools view           | N  |
| read 1 aligned by cycle         |                         | CIGAR                   | Y  |
| read 1 average length           |                         | samtools stats: FRL     | N  |
| read 1 deletion by cycle        |                         | CIGAR                   | Y  |
| read 1 hard clip by cycle       |                         | CIGAR                   | Y  |
| read 1 insertion by cycle       |                         | CIGAR                   | Y  |
| read 1 length histogram         |                         | samtools stats: FRL     | N  |
| read 1 mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read 1 quality by cycle         |                         | samtools stats: FFQ     | N  |
| read 1 quality histogram        |                         | samtools stats: FFQ     | N  |
| read 1 soft clip by cycle       |                         | CIGAR                   | Y  |
| read 2 aligned by cycle         |                         | CIGAR                   | Y  |
| read 2 average length           |                         | samtools stats: LRL     | N  |
| read 2 deletion by cycle        |                         | CIGAR                   | Y  |
| read 2 hard clip by cycle       |                         | CIGAR                   | Y  |
| read 2 insertion by cycle       |                         | CIGAR                   | Y  |
| read 2 length histogram         |                         | samtools stats: LRL     | N  |
| read 2 mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read 2 quality by cycle         |                         | samtools stats: LFQ     | N  |
| read 2 quality histogram        |                         | samtools stats: LFQ     | N  |
| read 2 soft clip by cycle       |                         | CIGAR                   | Y  |
| read ? aligned by cycle         |                         | CIGAR                   | Y  |
| read ? average length           |                         | CIGAR                   | Y  |
| read ? deletion by cycle        |                         | CIGAR                   | Y  |
| read ? hard clip by cycle       |                         | CIGAR                   | Y  |
| read ? insertion by cycle       |                         | CIGAR                   | Y  |
| read ? length histogram         |                         | CIGAR                   | Y  |
| read ? mismatch by cycle        |                         | samtools stats: MPC     | N  |
| read ? quality by cycle         |                         | CIGAR                   | Y  |
| read ? quality histogram        |                         | CIGAR                   | Y  |
| read ? soft clip by cycle       |                         | CIGAR                   | Y  |
| reads mapped and paired         | mate unmapped reads     | samtools stats: SN      | N  |
| reads on target                 |                         | bedtools                | Y  |
| reads per start point           |                         | RNAME/POS fields        | Y  |
| readsMissingMDtags              |                         | CIGAR                   | Y  |
| run name                        |                         | metadata                | N  |
| sample                          |                         | metadata                | N  |
| sample rate                     |                         | input parameter         | N  |
| soft clip bases                 |                         | CIGAR                   | Y  |
| target file                     |                         | input parameter         | N  |
| total reads                     |                         | samtools stats: SN      | N  |
| total target size               | target size             | bedtools                | N  |
| unmapped reads                  |                         | samtools stats: SN      | N  |

- The DS column indicates which metrics are affected by downsampling (if any).

## General notes on data sources

- `samtools` and `bedtools` are wrapped by the Python packages `pysam` and `pybedtools`, respectively.
- `pysam` is also used for processing CIGAR strings.
- `pairsMappedAbnormallyFar` is a secondary metric, computed from the insert size histogram, which in turn is derived from `samtools stats`. In this instance, 'abnormal' is defined by the `insert max` parameter.

## Filtering

There are two possible filter mechanisms in BamQC: _downsampling_ and _quality filtering_.

If both are in effect, downsampling is applied first.

### Downsampling

Downsampling is applied for faster and more efficient data processing. Large BAM files, with more than approximately 1 million reads, are downsampled at a rate of 1 in 1000. The downsampling rate is recorded in the JSON output as `sample rate`. 

The exact number of reads in a BAM file is not known unless one reads the entire file. To speed up decision-making, the number of reads is estimated from file size. A threshold of 100 MB corresponds to approximately 1 million reads of 2x151 length. Further details are on the [OICR wiki](https://wiki.oicr.on.ca/display/GSI/2019-04-04+Estimate+of+BAM+reads+from+file+size).

If downsampling is in effect, metrics derived from CIGAR strings and bedtools are evaluated on the downsampled reads. The estimated value of the metric on the full dataset must be scaled accordingly. For example, at a downsampling rate of 1 in 1000, if the reported 'hard clip bases' is 250, then there will be approximately 250,000 hard clip bases in the entire BAM file. Metrics derived directly from samtools will always be calculated on the full dataset, without downsampling. The `DS` column in the summary table records which metrics are affected by downsampling.

Downsampling is carried out using `samtools view` with a fixed random seed, so results will be consistent on repeated runs of bam-qc-metrics. For paired-end data, if a given read 1 is part of the downsampled set, its corresponding read 2 will be as well.

### Quality filtering

BamQC has an option to exclude reads with alignment quality below a given threshold.

If quality filtering is in effect, the `unmapped reads` metric is computed only on reads which _failed_ the filter -- because by definition, an unmapped read does not have an alignment quality. All other metrics are computed only on reads which _passed_ the filter.

## Special case: Sequence mismatches

The `mismatched bases` field may not be consistent with the mismatch-by-cycle fields for each read.

The former is derived from `samtools stats SN`; the latter from `samtools stats MPC` with a given alignment reference.
