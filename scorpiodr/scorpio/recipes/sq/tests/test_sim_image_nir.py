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
                 drive_url + "&id=1wofSNGArQbmaItnUKNOyIKsuHqbwJOBn",
                 "SCORPIO-H-2-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=1wMQB2E_9K8lDLNE9S3P7GcNcq8xXUL86",
                 "SCORPIO-H-3-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=13pgVAjconiw0LyIRyFu-fYQCSk4kXoxQ",
                 "SCORPIO-H-4-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=11QZBYehU2uZ6ZKR3iWRV4VeHgqvCWfs_",
                 "SCORPIO-H-5-IMAGING-FULL-DARK-20231211-142220-1.fits" :
                 drive_url + "&id=1REruNM1lgwakHLquFHDEQASSkaxufiAp"},
        "flat": {"SCORPIO-H-1-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1gz1M6xLc1kwhdWttN_Ot9tZRxxrHQKc3",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1xIBfS7VdzLMRQvT-SMXTJHn5hcCN_fW8",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1VRzOsi7AnT36hEOkl9kTyeZXTzZNzbrO",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1c6K3PpQVc3k_9-i7x28WwHa6iCzzBPjk",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DARK-20231121-094302-1.fits" :
                 drive_url + "&id=1IAJdHYXOE_UV_yzULrEwM9q0BtD9FEhJ",
                 "SCORPIO-H-1-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1Fh6E7TSKiYqrP55Jjr3TE1zh_h_4kQTU",
                 "SCORPIO-H-2-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1BrL94LXduxHpNsNIyDrLT9k41tNgjx0A",
                 "SCORPIO-H-3-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1rbuhfjetoAw6mat4SV_0VygcuNNNm0FZ",
                 "SCORPIO-H-4-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=1D4U4MEnGd2hKc4r7FelFO3szmmRvx2Hc",
                 "SCORPIO-H-5-IMAGING-FULL-FLAT-DOME-20231121-094038-1.fits" :
                 drive_url + "&id=15u9xgfuuCd-vg5maAEg-H1awbYT9zjvm"},
        "sci": {"SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-1.fits" :
                drive_url + "&id=1YCzT87JhOfPHCuyPK6K7qSJn9NpSvTAv",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-2.fits" :
                drive_url + "&id=1LHYfSt2TzYN9mXtuyqw8aqODmW3wMBHZ",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230-3.fits" :
                drive_url + "&id=1Bb0Btb8r3N3Fpc4PHwp-E8BB_8PBooEs",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-1.fits" :
                drive_url + "&id=10qq_PsiQqYcCv9KSaB0utnI4XXFpjO5S",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-2.fits" :
                drive_url + "&id=1iQ791_lPqOylbXUyyMXn7QfVDLhxabes",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129-3.fits" :
                drive_url + "&id=1FKyx0KNPl8AdP6MHloUUMrMgHwxrZ8NE",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-1.fits" :
                drive_url + "&id=1fMeTisDgtGj1Cl_9i8SAqJN-HfQ6z5Oc",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-2.fits" :
                drive_url + "&id=1Q02ylGVYGTgX-4mY6K0JrNuDYSqRxnS_",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155-3.fits" :
                drive_url + "&id=1l-E04G7JVFSWqFbX97nBCDD2TM2IDTHQ",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-1.fits" :
                drive_url + "&id=1xt02ttO0yzvuMJb4F8kcrxrHe6A1VG8q",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-2.fits" :
                drive_url + "&id=1JUecEKSpFD0OnYidVY6-p-Tt2wuhkjne",
                "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559-3.fits" :
                drive_url + "&id=1eyiklQioElJLHXE56iqoNz7_jBwmxq8Q"},
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
