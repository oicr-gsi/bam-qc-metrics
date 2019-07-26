CHANGELOG
=========

v0.1.3: 2019-07-26
------------------

Fixed:
- Use correct read length if all data is 'unknown read'
- Correctly handle empty read length histogram
- Correctly handle missing FFQ/LFQ in samtools stats

Added:
- 'Reads per start point' metric
- Sanity checks on output variables in Python tests

Changed:
- Run `samtools stats` _before_ downsampling and _after_ quality filtering (if any)
- Do downsampling using `samtools view` and a random seed instead of iterating over the reads
- Rename `trim_quality` parameter as `skip_below_mapq`
- Replace obsolete `distutils` with `setuptools` in setup.py


v0.1.2 : 2019-07-23
-------------------

Fixed:
- Count unmapped reads before, not after mapping quality filter is applied
- If read is aligned to reverse strand, also reverse order of CIGAR tuple
- Fixups for reading `mark_duplicates` text file

Added:
- Metric specifications document `metrics.md`

Changed:
- Initialize all instance variables at start of `__init__()` in `bam_qc.py`

v0.1.1 : 2019-07-12
-------------------

Fixed:
- Do not crash on empty CIGAR string

Added:
- Add install_requires in setup.py, allows pip to install prerequisites

Changed:
- Move `diff.py` script to `bin/` and rename as `diff_bam_qc_output.py`

v0.1.0 : 2019-07-11
-------------------

Initial release for testing

Added:
- Rewrite in Python of existing Perl code from https://github.com/oicr-gsi/bamqc
- Tests, setup.py and requirements.txt files
