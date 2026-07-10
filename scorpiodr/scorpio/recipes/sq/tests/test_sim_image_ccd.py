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
                 drive_url + "&id=1TV-fQgGC595AZhA8pzngLocMrULdzbZ0",
                 "SCORPIO-i-2-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=11OX11TGCxFBwfBKPP5vnlohDs3YUFYoF",
                 "SCORPIO-i-3-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1iW5VWzk-N4PZvyafSnG9Lc1wsb6IQeEK",
                 "SCORPIO-i-4-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1jU2Qv5J2fTpcseeZ1N8ZfHqPWiHoHmEo",
                 "SCORPIO-i-5-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1EyYFU499Bg6mGlovn9L4pNC45KiuKA-H"},
        "dark": {"SCORPIO-i-1-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1WM8Cgvm0Xp8lju_G1_HWKS2n6hcBKKAN",
                 "SCORPIO-i-2-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1IeZwy6azxcS6lfwYCYwgJZyCxDqDgsYP",
                 "SCORPIO-i-3-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1sqmug1PQC1WycxQuokZAMOCw85PppIde",
                 "SCORPIO-i-4-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1nL3GP5xqMTCLJ-qq31CYrJMizJU21OTr",
                 "SCORPIO-i-5-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1LPkqt1QvSSXUXO9t9XvZcdAw5v0VGuxb"},
        "flat": {"SCORPIO-i-1-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1StuLxJaCFF9XgBMhPMC1I6ae92tWpty6",
                 "SCORPIO-i-2-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1nuwlmmWHVrrO2oLy0UQXMMHWPUdfy4H8",
                 "SCORPIO-i-3-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1zuXBTE0gZplKNvTERQxzS0-QEbRes081",
                 "SCORPIO-i-4-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1FSTaQe2PvkDClDGEd8tIBmK0tKHM5Ej-",
                 "SCORPIO-i-5-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1FgKabAPgmKzoDEl4QRiQ2ZJM70GG7Gpe"},
        "sci": {"SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-1.fits" :
                drive_url + "&id=10PEUJ4wuGiODA5X7QwjYMnOwcQSMkuha",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-2.fits" :
                drive_url + "&id=1fvqcYYZLPknB8hlThjLanRWHFd_4ps2E",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-3.fits" :
                drive_url + "&id=1dUDCc0kQXYC6q1Wf1OVkNGJEr3g8oGlF",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-1.fits" :
                drive_url + "&id=1Bkp2soWvXRi5m17VWAzo91tKWyMCoXcS",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-2.fits" :
                drive_url + "&id=1k5FmtO2qDzbHI3eBBLHxRE-A961qye6O",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-3.fits" :
                drive_url + "&id=1ocnPDdLAT36mXihipTBQ7waI99SR07km",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-1.fits" :
                drive_url + "&id=1bFssqFy0gZK_7au7lU0H6GG6UAusAdGd",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-2.fits" :
                drive_url + "&id=1UTXqpruQaaAub4uMy9MLsAsxza2kXX0V",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-3.fits" :
                drive_url + "&id=1VjcEj3TPPmiXSGxHtAl3lRORz2kVGdZ0",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-1.fits" :
                drive_url + "&id=1OJeu4mD3LDv-aRxjgSAeNU0_jOX8dSie",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-2.fits" :
                drive_url + "&id=1p52d4L79SVvUFyIcvd86MoIvFBKWb4UB",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-3.fits" :
                drive_url + "&id=1I--GdUjDebS1TeOL5VmxKs84fpTxet3C"},
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
