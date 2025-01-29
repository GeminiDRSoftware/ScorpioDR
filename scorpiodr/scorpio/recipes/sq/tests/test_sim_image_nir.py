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
    "sim-H-imaging-full": {
        "dark": {"SCORPIO-H-1-IMAGING-FULL-DARK-20231211-142220.fits" :
                 drive_url + "&id=1I6eO8hrVBE3A3vaxP4HlrHeykiolOXfn",
                 "SCORPIO-H-2-IMAGING-FULL-DARK-20231211-142220.fits" :
                 drive_url + "&id=1u8mt6xHVJOLNDdeNlk6AzZAchk1ci2Md",
                 "SCORPIO-H-3-IMAGING-FULL-DARK-20231211-142220.fits" :
                 drive_url + "&id=1Yk2vxGX8KzYPL2TSfJdr1XZyjNlERR1C",
                 "SCORPIO-H-4-IMAGING-FULL-DARK-20231211-142220.fits" :
                 drive_url + "&id=14Evt-bMJCxvAZ0QXOwFDwORLcUm95a5z",
                 "SCORPIO-H-5-IMAGING-FULL-DARK-20231211-142220.fits" :
                 drive_url + "&id=1b5KM_vUZHOUUmWNUSTSP6RasZk69vpu5"},
        "flat": {"SCORPIO-H-1-IMAGING-FULL-FLAT-DARK-20231121-094302.fits" :
                 drive_url + "&id=1qA3cXg5Jc99ii6v84196CphrGHVLMa2F",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DARK-20231121-094302.fits" :
                 drive_url + "&id=1sJ70Gtl_xeD3NHw6Q_xJASYAIw5w-1R7",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DARK-20231121-094302.fits" :
                 drive_url + "&id=1n_p4b_qPoIS-jys4WxWCfddqq9Jc7srP",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DARK-20231121-094302.fits" :
                 drive_url + "&id=1iaF9VAsgXmI5np3mMEK0GNrLJnr3FDg4",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DARK-20231121-094302.fits" :
                 drive_url + "&id=1al2eOuuDWiXBHYWs2lOQccF9mky49NAa",
                 "SCORPIO-H-1-IMAGING-FULL-FLAT-DOME-20231121-094038.fits" :
                 drive_url + "&id=1Z645Js9yoVuOFG47l8_fNX1s2nP2KDof",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DOME-20231121-094038.fits" :
                 drive_url + "&id=1-ZiL1n1GcDpwXWLRo3aD5gAhk0K7ocQw",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DOME-20231121-094038.fits" :
                 drive_url + "&id=1Iz-CpylhDn8mRP3DsZ_8S7TGq5-3oYrc",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DOME-20231121-094038.fits" :
                 drive_url + "&id=1zUCyukK3JIfZkdY8yD-XDe8xARVG-81H",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DOME-20231121-094038.fits" :
                 drive_url + "&id=1TAyV7A12SdbkXc7E3euNzLt2mtSnHKnj"},
        "sci": {"SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230.fits" :
                drive_url + "&id=1DjNPm-VWWFXsltqJNOqTd9RlaZW8qbqL",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129.fits" :
                drive_url + "&id=1u-BugENXfOulu4nW6Wi1-dBDE0FXPwEE",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155.fits" :
                drive_url + "&id=1qp5FlPZyvGjiwbtlpVgR-2tGa5TV1rt3",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559.fits" :
                drive_url + "&id=1WMDYwvg8_XWqMe16VeE8I21T_JouI845"},
        "ucals": [],
        "refs" : ["SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230_stack.fits"],
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
        Keep pre-stack data? (Uses a lot of disk space; this currently doesn't
        apply for the NIR recipe, but the option is retained in case it's
        needed later).
    test_case : str
        Test parameter containing the key to the `datasets` dictionary.
    """
    with change_working_dir(test_case):

        cals = []

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
            run_reduce(sci_paths, f"sci_{test_case}", cals,
                       user_pars=datasets[test_case]["ucals"])
            # if not keep_data:
            #     print(' Deleting pre-stack files.')

        # Comparison of final output:
        for f in datasets[test_case]["refs"]:
            adout = adopen(f)
            adref = adopen(os.path.join(path_to_refs, f))
            assert ad_compare(adout, adref, atol=1e-4)


if __name__ == '__main__':
    pytest.main()
