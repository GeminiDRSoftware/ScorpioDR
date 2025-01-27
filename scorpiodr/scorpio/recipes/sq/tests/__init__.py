#!/usr/bin/env python
"""
Tests for SCORPIO SQ recipes.
"""

import pytest
import os

from recipe_system import __version__ as rs_version
from recipe_system.reduction.coreReduce import Reduce
from recipe_system.utils.reduce_utils import normalize_ucals, buildParser

from gempy.utils import logutils


# -- Helper functions (based on GMOS) -----------------------------------------
def run_reduce(file_list, label, calib_files, recipe_name=None, save_to=None,
               user_pars=None):
    """
    Helper function used to prevent replication of code.

    Parameters
    ----------
    file_list : list
        List of files that will be reduced.
    label : str
        Labed used on log files name.
    calib_files : list
        List of calibration files properly formatted for DRAGONS Reduce().
    recipe_name : str, optional
        Name of the recipe used to reduce the data.
    save_to : str, optional
        Stores the calibration files locally in a list.
    user_pars : list, optional
        List of user parameters

    Returns
    -------
    list : an updated list of calibration files
    """
    #objgraph = pytest.importorskip("objgraph")

    logutils.get_logger().info("\n\n\n")
    logutils.config(file_name=f"test_image_{label}.log")

    # Must instantiate with adpkg/drpkg up front to find the recipe:
    args = buildParser(rs_version).parse_args([])
    args.adpkg = "scorpio_instruments"  # above parse_args mishandles these!
    args.drpkg = "scorpiodr"            #
    args.files = file_list
    args.user_cal = calib_files
    args.userparam = user_pars
    r = Reduce(sys_args=args)

    if recipe_name:
        r.recipename = recipe_name

    r.runr()

    if save_to:
        calib_files.append("{}:{}".format(
            save_to, os.path.join("calibrations", save_to, r._output_filenames[0])))
        [os.remove(f) for f in r._output_filenames]

    # check that we are not leaking objects
    #assert len(objgraph.by_type('NDAstroData')) == 0

    return calib_files


# -- Custom configuration -----------------------------------------------------
@pytest.fixture(scope='module')
def keep_data(request):
    """
    By default, the tests will delete pre-stack files to save disk space. If one
    needs to keep them for debugging, one can pass --keep-data argument to the
    command line call to force the tests to keep this data.

    Parameters
    ----------
    request : fixture
        Represents the test that calls this function.
    """
    return request.config.getoption("--keep-data")
