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
                 drive_url + "&id=1taeCegHnI5IGJGoqM1xqbQlvJ04ebKoQ",
                 "SCORPIO-i-2-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1bYOYLne_eDFOPC2RQSDkqqBVKT2ZmHuo",
                 "SCORPIO-i-3-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=10grBEPBUY-IFzvYM88WIxNnrNPAtL72O",
                 "SCORPIO-i-4-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1p7u7lEibkffjdb938hDCBqyFzkJrOdK7",
                 "SCORPIO-i-5-IMAGING-FULL-BIAS-20231003-130052-1.fits" :
                 drive_url + "&id=1uDljQxtrmC5s7g6ECqCipeSEKaUPC2w9"},
        "dark": {"SCORPIO-i-1-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1iSXN17Gq2AKNihYUxo_yoIH4B5XcFt4h",
                 "SCORPIO-i-2-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1cUP46gAG_vefP8jppB_aMPhasDXPXKYj",
                 "SCORPIO-i-3-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1zxmx_bVDLzxMe9uNPIOtmwXagg9SU5eh",
                 "SCORPIO-i-4-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1TVS7aBs2ab65DAlgFPGtNEC464lmpBxS",
                 "SCORPIO-i-5-IMAGING-FULL-DARK-20231003-130122-1.fits" :
                 drive_url + "&id=1bYqZgYezK-BusnnJyvMm0lVrgecE0zM1"},
        "flat": {"SCORPIO-i-1-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1H1kQjgmnSmTytVRa3GQzRh0-r6nXWfhB",
                 "SCORPIO-i-2-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=13yyZyRtKNlp7MjUsxzFPThGeWmFSh5jX",
                 "SCORPIO-i-3-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1B8lzVogvqj_lf06ktxZsJTiA_De4igv9",
                 "SCORPIO-i-4-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1vM75c74PxL8_xe4pyc3bo4MZnNUG1p5K",
                 "SCORPIO-i-5-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits" :
                 drive_url + "&id=1FiD2dZjp5AbF3pSPe3JRj9jLIgrBXMvw"},
        "sci": {"SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-1.fits" :
                drive_url + "&id=1eINdo2Vj41Yyk99tWb9EfNJtwXRs2nWp",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-2.fits" :
                drive_url + "&id=1imxyNudBLc8VJeSi2irOkbNJCEDtCIZg",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-3.fits" :
                drive_url + "&id=1i_Ez4DN0go_r3RdaMXqPjG23oz_4Wf6P",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-1.fits" :
                drive_url + "&id=1NXxJEBtRS_FDWt1eGsSFY-vy9g6eY6gY",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-2.fits" :
                drive_url + "&id=1mV1vXurgDt6NGMretOvyHG_96iRikNBC",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-161922-3.fits" :
                drive_url + "&id=1VK9f3OKuzzISJEdbyXMgM_7O9A5M122B",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-1.fits" :
                drive_url + "&id=1GkqpB1WFU7jpxtJBaAteTfgnLrNRBXtf",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-2.fits" :
                drive_url + "&id=10AIRoU6lGXFhgOrJFdkMV8NwDnttfWnt",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-163737-3.fits" :
                drive_url + "&id=1WPCigwCZK3-JLoN2lxqRVIxPXYyoSqEW",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-1.fits" :
                drive_url + "&id=11wOAznHzJxjB02Wo1mffUF-PLepBbk5p",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-2.fits" :
                drive_url + "&id=1jwVtipssOVB2AVFl-GPUnJ-FzSPyhs90",
                "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240917-080911-3.fits" :
                drive_url + "&id=1pw83gbwXUJ9-5FfIOTC058GolghGyHCK"},
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
