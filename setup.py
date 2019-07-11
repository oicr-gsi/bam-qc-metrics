"""
Setup script for bam-qc-metrics
"""

import sys
from distutils.core import setup

if sys.version_info[0] != 3:
    print("ERROR: Python3 is required.", file=sys.stderr)
    sys.exit(1)

def main():
    setup(
        name = "bam-qc-metrics",
        version = "0.1.0",
        scripts = ['bin/run_bam_qc.py',],
        packages = ['bam_qc_metrics'],
        requires = ['pybedtools', 'pysam'],
        python_requires='>=3.5',
        author = "Iain Bancarz",
	author_email = "ibancarz@oicr.on.ca",
        description = "BAM QC metrics",
        long_description = "Python implementation of BAM QC metrics in use at the Ontario Institute for Cancer Research, https://oicr.on.ca",
        url = "https://github.com/oicr-gsi/bam-qc-metrics",
    )




if __name__ == "__main__":
	main()
