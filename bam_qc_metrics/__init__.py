
from .bam_qc import bam_qc, fast_metric_finder

# fast_metric_finder class is exported for tests

import os

def read_package_version():
    # This method depends on relative path to the VERSION file
    # So it has been placed in the package __init__ file, whose location will never change
    in_path = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir, 'VERSION'))
    with open(in_path) as version_file:
        package_version = version_file.read().strip()
    return package_version
