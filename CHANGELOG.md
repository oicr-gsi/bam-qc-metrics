CHANGELOG
=========

v0.1.1 : 2019-07-12
-------------------

- Do not crash on empty CIGAR string
- Add install_requires in setup.py, allows pip to install prerequisites
- Remove requirements.txt, unnecessary given improved setup.py
- Move `diff.py` script to `bin/` and rename as `diff_bam_qc_output.py`

v0.1.0 : 2019-07-11
-------------------

Initial release for testing

Added:
- Rewrite in Python of existing Perl code from https://github.com/oicr-gsi/bamqc
- Tests, setup.py and requirements.txt files
