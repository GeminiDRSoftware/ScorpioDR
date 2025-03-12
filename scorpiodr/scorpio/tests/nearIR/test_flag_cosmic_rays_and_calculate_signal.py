import os
import pathlib
import sys

import numpy as np
import pytest

from astrodata import open as adopen
from astrodata.testing import ad_compare, download_from_archive
from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.utils import logutils
from recipe_system.testing import ref_ad_factory

import scorpio_instruments
from scorpiodr.scorpio.primitives_scorpio_nearIR import ScorpioNearIR

from scorpiodr.scorpio.tests import ad

datasets = [
    "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230_refpixelsCorrected.fits",
#   "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129_refpixelsCorrected.fits",
#   "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155_refpixelsCorrected.fits",
    "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559_refpixelsCorrected.fits",
]

@pytest.mark.scorpioimage
@pytest.mark.preprocessed_data
@pytest.mark.regression
@pytest.mark.parametrize("ad", datasets, indirect=True)
def test_regression_in_flag_cosmic_rays_and_calculate_signal(
        ad, change_working_dir, ref_ad_factory
):
    with change_working_dir():
        logutils.config(
            file_name='log_regression_{:s}.txt'.format(ad.data_label())
        )
        p = ScorpioNearIR([ad])
        p.flagCosmicRaysFromNDRs()
        p.calculateSignalByRegression()
        out_ad = p.writeOutputs().pop()

    ref_ad = ref_ad_factory(out_ad.filename)

    for out_ext, ref_ext in zip(out_ad, ref_ad):
        np.testing.assert_allclose(out_ext.data, ref_ext.data, atol=1e-4)
        np.testing.assert_equal(out_ext.mask & DQ.cosmic_ray,
                                ref_ext.mask & DQ.cosmic_ray)


# -- Recipe to create pre-processed data --------------------------------------
def create_inputs_recipe():
    """
    Creates input data for tests.

    The raw files will be downloaded and saved inside the path stored in the
    `$DRAGONS_TEST/raw_inputs` directory. Processed files will be stored inside
    a new folder called "dragons_test_inputs". The sub-directory structure
    should reflect the one returned by the `path_to_inputs` fixture.
    """

    drive_url = "https://docs.google.com/uc?export=download"

    fnames = {
        "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-103230.fits" :
        drive_url + "&id=1DjNPm-VWWFXsltqJNOqTd9RlaZW8qbqL",
#       "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-115129.fits" :
#       drive_url + "&id=1u-BugENXfOulu4nW6Wi1-dBDE0FXPwEE",
#       "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240918-173155.fits" :
#       drive_url + "&id=1qp5FlPZyvGjiwbtlpVgR-2tGa5TV1rt3",
        "SCORPIO-H-1-IMAGING-FULL-OBJECT-20240919-094559.fits" :
        drive_url + "&id=1WMDYwvg8_XWqMe16VeE8I21T_JouI845",
    }

    path = pathlib.Path('dragons_test_inputs')
    path = (path / "scorpiodr" / "scorpio" / "nearIR" /
            "test_flag_cosmic_rays_and_calculate_signal" / "inputs")
    path.mkdir(exist_ok=True, parents=True)
    os.chdir(path)
    print('Current working directory:\n    {!s}'.format(path.cwd()))

    for fname, url in fnames.items():
        sci_ad = adopen(download_from_archive(fname, url=url))
        data_label = sci_ad.data_label()

        print('===== Reducing pre-processed data =====')
        logutils.config(file_name=f'log_{data_label}.txt')
        p = ScorpioNearIR([sci_ad])
        p.prepare()
        p.addDQ()
        p.ADUToElectrons()
        p.addVAR(read_noise=True, poisson_noise=True)
        p.nonlinearityCorrect()
        p.referencePixelsCorrect()
        p.writeOutputs()


if __name__ == '__main__':
    if "--create-inputs" in sys.argv[1:]:
        create_inputs_recipe()
    else:
        pytest.main([__file__])

