#!/usr/bin/env python

import glob
import pytest
import os

from astrodata import open as adopen
from astrodata.testing import ad_compare, download_from_archive
from recipe_system import __version__ as rs_version
from recipe_system.reduction.coreReduce import Reduce
from recipe_system.utils.reduce_utils import normalize_ucals, buildParser

from gempy.utils import logutils

drive_url = "https://docs.google.com/uc?export=download"

datasets = {

    # Simulated dataset doesn't follow the usual DHS naming:
    "sim-i-imaging-full": {
        "bias": {"SCORPIO-i-1-IMAGING-FULL-BIAS-20231003-130052.fits" :
                 drive_url + "&id=1Qc69BG_nq1U3f1a7k36804BUwcN8-yip",
                 "SCORPIO-i-2-IMAGING-FULL-BIAS-20231003-130052.fits" :
                 drive_url + "&id=1ipJ4zWM-BecD6ef1iyB7ejx-b3a_jMZ2",
                 "SCORPIO-i-3-IMAGING-FULL-BIAS-20231003-130052.fits" :
                 drive_url + "&id=1rMnEGwkuMAqCYJ2UkZV4cqIDSxYppWwJ",
                 "SCORPIO-i-4-IMAGING-FULL-BIAS-20231003-130052.fits" :
                 drive_url + "&id=1mbZWp7EJkvWitkSMeQ-O_NwIJJ5EPGHK",
                 "SCORPIO-i-5-IMAGING-FULL-BIAS-20231003-130052.fits" :
                 drive_url + "&id=1m4IRGjRSDBWIYP2Yt8IqYnq4v7xY0BO2"},
        "dark": {"SCORPIO-i-1-IMAGING-FULL-DARK-20231003-130122.fits" :
                 drive_url + "&id=1ulJdpJE20i6Y4QPNPoW-GW36Vx_S9eJp",
                 "SCORPIO-i-2-IMAGING-FULL-DARK-20231003-130122.fits" :
                 drive_url + "&id=15A4Qc7I8xBFj-7qPzclEzgXGkC5ZKcRQ",
                 "SCORPIO-i-3-IMAGING-FULL-DARK-20231003-130122.fits" :
                 drive_url + "&id=1v9GkdKoKqABQlcgqmX6vXXZSmP8JOIMJ",
                 "SCORPIO-i-4-IMAGING-FULL-DARK-20231003-130122.fits" :
                 drive_url + "&id=1SNN4eJ0YInwUfCedaFq_QaSswLhqiyuS",
                 "SCORPIO-i-5-IMAGING-FULL-DARK-20231003-130122.fits" :
                 drive_url + "&id=1nqzTufnz4y0ch0Bhnc-q0812ZimPzeVm"},
        "flat": {"SCORPIO-i-1-IMAGING-FULL-FLAT-DOME-20231211-120210.fits" :
                 drive_url + "&id=1C89vnsEczFWki8hHv3AifXu1Um-WMC3d",
                 "SCORPIO-i-2-IMAGING-FULL-FLAT-DOME-20231211-120210.fits" :
                 drive_url + "&id=1fluTIFvbTwm3yeZ_TAlmRfRp4HrIfmo5",
                 "SCORPIO-i-3-IMAGING-FULL-FLAT-DOME-20231211-120210.fits" :
                 drive_url + "&id=1En3iW8o65LeVz637yatmaaQKDq0GZAff",
                 "SCORPIO-i-4-IMAGING-FULL-FLAT-DOME-20231211-120210.fits" :
                 drive_url + "&id=1GzXK7msWiuvJ7HcjpIYRJ50GYE0whQxJ",
                 "SCORPIO-i-5-IMAGING-FULL-FLAT-DOME-20231211-120210.fits" :
                 drive_url + "&id=1DYbyBIGijaaLRzlMadPscFDHECcgn-mK"},
        "sci": {"SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030.fits" :
                drive_url + "&id=1pLRgLyptTsIAeQEjS_xz2Wjw8CtSJcAo",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922.fits" :
                drive_url + "&id=1yjs026Kl4O5ckVLO4MgXfcIuqoftOq4e",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737.fits" :
                drive_url + "&id=1AS1r2J_kiENTYEU21mU0-hRogJZCzvZ6",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911.fits" :
                drive_url + "&id=1AiNhzVYGbp3fcw064t09b3HEjCIOleGt"},
        "ucals": [],
    },

}

@pytest.mark.scorpioim
@pytest.mark.integration_test
@pytest.mark.dragons_remote_data
@pytest.mark.parametrize("test_case", list(datasets.keys())[:1])
def test_reduce_image(change_working_dir, keep_data, test_case, path_to_refs):
    """
    Tests that we can run all the data reduction steps on a complete dataset.

    Parameters
    ----------
    change_working_dir : fixture
        Change the current work directory using a context manager.
    keep_data : fixture
        Keep pre-stack data? (Uses a lot of disk space)
    test_case : str
        Test parameter containing the key to the `datasets` dictionary.
    """
    with change_working_dir(test_case):

        cals = []

        # Biases
        bias_paths = [download_from_archive(f, url=u) for
                      f, u in datasets[test_case]["bias"].items()]
        cals = run_reduce(bias_paths, f"bias_{test_case}", cals,
                          save_to="processed_bias")

        # Darks
        bias_paths = [download_from_archive(f, url=u) for
                      f, u in datasets[test_case]["dark"].items()]
        cals = run_reduce(bias_paths, f"bias_{test_case}", cals,
                          save_to="processed_dark")

        # Flats
        flat_paths = [download_from_archive(f, url=u) for
                      f, u in datasets[test_case]["flat"].items()]
        cals = run_reduce(flat_paths, f"flat_{test_case}", cals,
                          save_to="processed_flat")

        # Science data
        if "sci" in datasets[test_case]:
            sci_paths = [download_from_archive(f, url=u) for
                         f, u in datasets[test_case]["sci"].items()]
            # cals = run_reduce(
            #     sci_paths, f"fringe_{test_case}", cals,
            #     recipe_name='makeProcessedFringe',
            #     save_to="processed_fringe")
            run_reduce(sci_paths, f"sci_{test_case}", cals,
                       user_pars=datasets[test_case]["ucals"])
            if not keep_data:
                print(' Deleting pre-stack files.')
                [os.remove(f) for f in glob.glob("*_CRMasked.fits")]
                [os.remove(f) for f in glob.glob("*_align.fits")]

        # Comparison of final output:
        for f in ['SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030_stack.fits']:
            adout = adopen(f)
            adref = adopen(os.path.join(path_to_refs, f))
            assert ad_compare(adref, adout, atol=1e-4)


# -- Helper functions ---------------------------------------------------------
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


if __name__ == '__main__':
    pytest.main()
