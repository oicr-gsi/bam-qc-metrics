CHANGELOG
=========

Unreleased
----------

Fixed:
- GP-2055 Run without error if `--target` not given

Added:
- GP-2045 Implementation of additional bedtools metrics
- GP-2048 Downsampling to exact number of reads
- Command-line options `--all-reads` and `--sample` control sample level
- Default number of reads to sample is 1 million
- Downsampling has no effect if fewer than 1 million reads present
- GP-2057 Logging and profiling
- Command-line options for `--log` and `--profile`
- `--verbose` and `--debug` options set the log level
- Option to set random seed, with updated tests
- GP-2093 Bedtools coverage depth metrics
- Additional metrics: Total coverage, coverage per target, target size, coverage histogram

v0.1.6: 2019-08-16
------------------

Fixed:
- Bugfix: `setup.py` now copies the VERSION file to the installation directory

v0.1.5: 2019-08-16
------------------

Fixed:
- GP-2035 More efficient finding unmapped reads
- GP-2036 Failure to skip BED file headers
- GP-2037 Count mismatches using MD tag, not CIGAR string
- GP-2042 Correct reads per start point calculation

Added:
- GP-2038 Report software versions in JSON
- Package version is read from the VERSION file
- Specify workflow version using the '-w' option for `run_bam_qc.py`
- Convenience script `update_test_data_version.py` to update package version in test data
- GP-2045 Placeholders for additional bedtools metrics

Changed:
- GP-2040 refactoring
- Separate classes for "fast" metrics (found before downsampling) and "slow" (after)

v0.1.4: 2019-07-26
------------------

Fixed:
- Bug in downsampling; `samtools view` expects downsample rate as decimal, not integer


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
