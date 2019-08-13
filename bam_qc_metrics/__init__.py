
from .bam_qc import bam_qc, fast_metric_finder, version_updater

# 'base' and 'fast_metric_finder' classes are exported for tests
# 'version_updater' class is exported for update_test_data_version.py

import os

def read_package_version():
    # This method depends on relative path to the VERSION file
    # So it has been placed in the package __init__ file, whose location will never change
    in_path = os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir, 'VERSION'))
    with open(in_path) as version_file:
        package_version = version_file.read().strip()
    return package_version

def get_data_dir_path():
    # return the path of the 'data' directory
    return os.path.realpath(os.path.join(os.path.dirname(__file__), os.pardir, 'data'))
