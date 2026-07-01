#!/usr/bin/env python

import glob
import pytest
import os

from astrodata import open as adopen
from astrodata.testing import ad_compare, download_from_archive

from . import run_reduce, keep_data

drive_url = "https://docs.google.com/uc?export=download"

datasets = {

    # Simulated dataset doesn't follow the usual DHS naming:
    "sim-i-imaging-full": {
        "bias": {"SCORPIO-i-1-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1reUybEi-4bpfCvA-TVTDL2I9OoI-K7xN",
                 "SCORPIO-i-2-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1GvF5rquFbIC1b9uiRvUCKt4mHgAiyj5v",
                 "SCORPIO-i-3-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1kLK57n5PbT00GcmyGa5snj2NiFcarg3H",
                 "SCORPIO-i-4-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1FGhKyvt5f5sAHyJweT8pQ0hGrOFWxERU",
                 "SCORPIO-i-5-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1UapbeT3YcjC4Iw8Avsrk3c0w6ijmv2eq"},
        "dark": {"SCORPIO-i-1-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1c1ANFZg-cqqDAGxNyFP-_cYEybFyHbLg",
                 "SCORPIO-i-2-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1dnBFECKvanNQebZIlJr6LkE5AaUXmb8e",
                 "SCORPIO-i-3-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1FH03OeUVyLWVrWA4TKUqkpCU4G9l8ot3",
                 "SCORPIO-i-4-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1fMQzDeZNUE0Jd5wlNVf53mrEOcupc8wC",
                 "SCORPIO-i-5-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=18eCnyrJyiXb1LMT_MB47hYInZMJhzu1H"},
        "flat": {"SCORPIO-i-1-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1JP43k_fNXndUSZOSti2v3aD78xguL4_P",
                 "SCORPIO-i-2-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=13yivmWNll3189-2pMI3EaPBBhmFYduih",
                 "SCORPIO-i-3-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1vTVcFwCZASpjggfEWijr3efqh-Amz7GL",
                 "SCORPIO-i-4-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1ToDSE_CJS7k2l4KaStINt5XJADidumKZ",
                 "SCORPIO-i-5-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1J6PpJNsjG-jikj7EdLuMmb97plbsLQTc"},
        "sci": {"SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-1.fits" :
                drive_url + "&id=1x1EiEHP4L9VAnWJgxS7zkuxzt2q7BCBt",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-2.fits" :
                drive_url + "&id=1ZvlIr3k90iRKareoxNmSyN_onuqTfYh4",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-3.fits" :
                drive_url + "&id=1woyFepT8XNy-TBrJqVwAV3op9rtqqdmT",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-1.fits" :
                drive_url + "&id=1r2q3DMBf7f7tOsQQyjijBv8RXMM-0Woj",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-2.fits" :
                drive_url + "&id=18dXzWTPzVHnCQnTuUMXGctj8BpxNZFOU",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-3.fits" :
                drive_url + "&id=1GvF5rquFbIC1b9uiRvUCKt4mHgAiyj5v",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-1.fits" :
                drive_url + "&id=1omRnc8KY0Y9mWY1zxkl-nlQtJlItqADu",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-2.fits" :
                drive_url + "&id=100zudHTtnQeo6S-NinihK1ciwF1Up9Ju",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-3.fits" :
                drive_url + "&id=10lV6Dq_o4-6T41vm_0_zOPmVimCGdihI",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-1.fits" :
                drive_url + "&id=1FM_bLSotWZWUifGdomnzXKlxA0uEnxVn",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-2.fits" :
                drive_url + "&id=1cknSYQ4377-TUtXuaWvOJ0XjaO737ECq",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-3.fits" :
                drive_url + "&id=1JQwNO8UV2hHRn7qp2WK16uYzMgigrrsc"},
        "ucals": [],
        "refs" : ["SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-1_stack.fits"],
    },

}

@pytest.mark.scorpioimage
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
        for f in datasets[test_case]["refs"]:
            adout = adopen(f)
            adref = adopen(os.path.join(path_to_refs, f))
            assert ad_compare(adout, adref, atol=1e-4)


if __name__ == '__main__':
    pytest.main()
