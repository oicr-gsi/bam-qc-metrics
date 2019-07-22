# Metric specifications for BamQC

BamQC metrics are computed using a number of methods, including third-party software tools, and output in JSON format. The output also contains metadata, such as the instrument and lane names; and input parameters, such as the target file. This document summarizes how metrics are computed, and details some special cases of interest.

Some metrics have been renamed in the move from the old Perl to new Python implementation.

## Summary of fields

| Name                            | Former name             | Source                  |
| --------------------------------|-------------------------|-------------------------|
| average read length             |                         | samtools stats: RL      |
| barcode                         |                         | metadata                |
| bases mapped                    | aligned bases           | samtools stats: SN      |
| deleted bases                   |                         | samtools stats: ID      |
| hard clip bases                 |                         | CIGAR                   |
| insert max                      |                         | input parameter         |
| insert size average             | insert mean             | samtools stats: SN      |
| insert size histogram           | insert histogram        | samtools stats: IS      |
| insert size standard deviation  | insert stdev            | samtools stats: SN      |
| inserted bases                  |                         | samtools stats: ID      |
| instrument                      |                         | metadata                |
| lane                            |                         | metadata                |
| library                         |                         | metadata                |
| mapped reads                    |                         | samtools stats: SN      |
| mark duplicates                 |                         | Picard MarkDuplicates   |
| mismatched bases                | mismatch bases          | samtools stats: SN      |
| non primary reads               |                         | samtools stats: SN      |
| number of targets               |                         | bedtools                |
| paired end                      | number of ends          | samtools stats: SN      |
| paired reads                    |                         | samtools stats: SN      |
| pairsMappedAbnormallyFar        |                         | insert size histogram   |
| pairsMappedToDifferentChr       |                         | samtools stats: SN      |
| properly paired reads           |                         | samtools stats: SN      |
| qual cut                        |                         | input parameter         |
| qual fail reads                 |                         | samtools view           |
| read 1 aligned by cycle         |                         | CIGAR                   |
| read 1 average length           |                         | samtools stats: FRL     |
| read 1 deletion by cycle        |                         | CIGAR                   |
| read 1 hard clip by cycle       |                         | CIGAR                   |
| read 1 insertion by cycle       |                         | CIGAR                   |
| read 1 length histogram         |                         | samtools stats: FRL     |
| read 1 mismatch by cycle        |                         | CIGAR                   |
| read 1 quality by cycle         |                         | samtools stats: FFQ     |
| read 1 quality histogram        |                         | samtools stats: FFQ     |
| read 1 soft clip by cycle       |                         | CIGAR                   |
| read 2 aligned by cycle         |                         | CIGAR                   |
| read 2 average length           |                         | samtools stats: LRL     |
| read 2 deletion by cycle        |                         | CIGAR                   |
| read 2 hard clip by cycle       |                         | CIGAR                   |
| read 2 insertion by cycle       |                         | CIGAR                   |
| read 2 length histogram         |                         | samtools stats: LRL     |
| read 2 mismatch by cycle        |                         | CIGAR                   |
| read 2 quality by cycle         |                         | samtools stats: LFQ     |
| read 2 quality histogram        |                         | samtools stats: LFQ     |
| read 2 soft clip by cycle       |                         | CIGAR                   |
| read ? aligned by cycle         |                         | CIGAR                   |
| read ? average length           |                         | CIGAR                   |
| read ? deletion by cycle        |                         | CIGAR                   |
| read ? hard clip by cycle       |                         | CIGAR                   |
| read ? insertion by cycle       |                         | CIGAR                   |
| read ? length histogram         |                         | CIGAR                   |
| read ? mismatch by cycle        |                         | CIGAR                   |
| read ? quality by cycle         |                         | CIGAR                   |
| read ? quality histogram        |                         | CIGAR                   |
| read ? soft clip by cycle       |                         | CIGAR                   |
| reads mapped and paired         | mate unmapped reads     | samtools stats: SN      |
| reads on target                 |                         | bedtools                |
| readsMissingMDtags              |                         | CIGAR                   |
| run name                        |                         | metadata                |
| sample                          |                         | metadata                |
| sample rate                     |                         | input parameter         |
| soft clip bases                 |                         | CIGAR                   |
| target file                     |                         | input parameter         |
| total reads                     |                         | samtools stats: SN      |
| total target size               | target size             | bedtools                |
| unmapped reads                  |                         | samtools stats: SN      |

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

If downsampling is in effect, all metrics are reported on the downsampled reads only. For example, if we have downsampled 1 in 1000 reads, the number of deleted bases will be approximately 1000 times fewer than in the full-sized BAM file.

Downsampling is deterministic; there is no random element. At a downsampling rate of 1000, every 2000th _pair_ of reads from the name-sorted BAM file is included. For paired-end data, this should ensure that the paired read 1 and read 2 are downsampled.

### Quality filtering

BamQC has an option to exclude reads with alignment quality below a given threshold.

If quality filtering is in effect, the `unmapped reads` metric is computed only on reads which _failed_ the filter -- because by definition, an unmapped read does not have an alignment quality. All other metrics are computed only on reads which _passed_ the filter.

## Special case: Sequence mismatches

The `mismatched bases` field may not be consistent with the mismatch-by-cycle fields for each read.

The former is derived from `samtools stats`; the latter from CIGAR strings. The CIGAR operation M for 'alignment match' may represent a sequence match or mismatch. So, the CIGAR string does not necessarily record all mismatches. In the event of inconsistency, `mismatched bases` should be taken as correct.
