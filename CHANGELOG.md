CHANGELOG
=========

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
