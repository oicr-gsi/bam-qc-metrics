"""
Setup script for bam-qc-metrics
"""

import os
from setuptools import setup

with open(os.path.join(os.path.dirname(__file__), 'VERSION')) as version_file:
    package_version = version_file.read().strip()

setup(
    name="bam-qc-metrics",
    version=package_version,
    scripts=['bin/run_bam_qc.py', ],
    packages=['bam_qc_metrics'],
    install_requires=['attrs', 'jsonschema', 'pybedtools', 'pyrsistent', 'pysam', 'six'],
    python_requires='>=3.5',
    author="Iain Bancarz",
    author_email="ibancarz@oicr.on.ca",
    description="BAM QC metrics",
    long_description="Python implementation of BAM QC metrics in use at the Ontario Institute for Cancer Research, https://oicr.on.ca",
    url="https://github.com/oicr-gsi/bam-qc-metrics",
)
