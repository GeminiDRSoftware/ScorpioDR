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
        "dark": {"SCORPIO-H-1-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=13eVw7Za38agX0pGukjxV-NQ6IteJBgY7",
                 "SCORPIO-H-2-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=1WKNOGGGP3OgXCPltOJPfR3dR__uDIi3i",
                 "SCORPIO-H-3-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=1qGiW3Wl7eKbHlQBxpWte4BawMzrRRzPD",
                 "SCORPIO-H-4-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=1__1r2D18US24rcHm_LzduBLwUotflM4f",
                 "SCORPIO-H-5-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=16TGZBaIKqeaFF4Th4NvtLw0nm6WoZiib"},
        "flat": {"SCORPIO-H-1-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1UBXK4PcMDSfDvM67HPKp8W83NrvQJEWr",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=15bHi2IKd109lrOs_mq19ghQy2dMFWLMj",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=14AlXEIVihxPubDo9Kk2qnHXsSs2HHiGL",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1l15dMejMNTeXXdPMqs6CIeWSgZHqhBNj",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=12AscAO0BsJqsrrp8KShjxxWOmMGvMzH_",
                 "SCORPIO-H-1-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1fn2FyYOPSI9p8P1EiUhQglgo_VTRnR6a",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1vOtSP6WEoB70eqLVlsZN_eOso7C6XNyh",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1gxG-2m8AJiB9af-EAYkaz1sTl0ZTaeUH",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1kNoHXhrKZelae1YSpuWlilrW8bxa8_RB",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=18hujKa1HfD4L48KBUvFadQJZ4kROIYO-"},
        "sci": {"SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-1.fits" :
                drive_url + "&id=1tHl3iKegDWWy0UDZrzV0pUZ5mH00Ekdm",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-2.fits" :
                drive_url + "&id=10WocJHl81eyfgOg6Zzh8hOmJ34NRqyyL",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-3.fits" :
                drive_url + "&id=1QpXPq14z7PMkRJVVZizFRMmNQkGeT9d2",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-1.fits" :
                drive_url + "&id=174d-851UcgWaHxBz-pKlhQQLQBVQatlI",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-2.fits" :
                drive_url + "&id=1O0U3-Ax2ZgPvsKYMK6edGFu16HD6uuDO",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-3.fits" :
                drive_url + "&id=1p76L8D9PZ8imMVtV1e03h3ohxwGaM-DV",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-1.fits" :
                drive_url + "&id=1LKrEuAEHZ-gCNpTWmAVOHawTXCK4-FI4",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-2.fits" :
                drive_url + "&id=1Bep_Z29OpYYIW45M812m9Lodo1KH0gti",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-3.fits" :
                drive_url + "&id=1GFBRxw5AHr8x9Pwc9BBWMTKZKLHYJSLV",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-1.fits" :
                drive_url + "&id=18oVmxvGrWVnh0SlryWcxmaRuZSY-iXa7",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-2.fits" :
                drive_url + "&id=1MEMhdNKH6cz_dHzZtW0y3YuxXhG-PtHh",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-3.fits" :
                drive_url + "&id=1nVEdINR-irCLGrbf9Yp3xH2plfvsZXSt"},
        "ucals": [],
        "refs" : ["SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-1_stack.fits"],
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
