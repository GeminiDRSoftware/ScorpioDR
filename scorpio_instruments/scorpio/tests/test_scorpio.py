import pytest

import astrodata
import astrodata.testing
import scorpio_instruments

import numpy as np


SCORPIO_DESCRIPTORS_TYPES = [
    ('detector_x_offset', float),
    ('detector_y_offset', float),
    ('pixel_scale', float),
  # ('nod_count', tuple),
  # ('nod_offsets', tuple),
  # ('shuffle_pixels', int),
]

test_files = [
    # Simulated images, pending real data:
    "SCORPIO-i-1-IMAGING-FULL-BIAS-20231003-130052-1.fits",
    "SCORPIO-i-1-IMAGING-FULL-DARK-20231003-130122-1.fits",
    "SCORPIO-i-1-IMAGING-FULL-FLAT-DOME-20231211-120210-1.fits",
    "SCORPIO-i-1-IMAGING-FULL-OBJECT-20240916-155030-1.fits",
]


@pytest.fixture(params=test_files)
def ad(request):
    filename = request.param
    path = astrodata.testing.download_from_archive(
        filename, url="https://docs.google.com/uc?export=download"
    )
    return astrodata.open(path)


@pytest.mark.dragons_remote_data
def test_is_right_instance(ad):
    assert isinstance(ad, scorpio_instruments.scorpio.adclass.AstroDataScorpio)


@pytest.mark.dragons_remote_data
def test_can_return_instrument(ad):
    assert ad.phu['INSTRUME'] == 'SCORPIO'
    assert ad.instrument() == ad.phu['INSTRUME']


@pytest.mark.dragons_remote_data
def test_can_return_ad_length(ad):
    assert len(ad)


@pytest.mark.parametrize("descriptor,expected_type", SCORPIO_DESCRIPTORS_TYPES)
@pytest.mark.dragons_remote_data
def test_descriptor_matches_type(ad, descriptor, expected_type):
    value = getattr(ad, descriptor)()
    assert isinstance(value, expected_type) or value is None, \
        "Assertion failed for file: {}".format(ad.filename)


# def test_tag_as_standard_fake(astrofaker):
#     # LTT4363 (a high proper motion specphot) on Jan 1, 2021
#     ad = astrofaker.create('SCORPIO', ['SPECT'],
#                            extra_keywords={'RA': 176.46534847,
#                                            'DEC': -64.84352513,
#                                            'DATE-OBS': '2021-01-01T12:00:00.000',
#                                            'OBSTYPE': 'OBJECT'}
#                            )
#     assert 'STANDARD' in ad.tags


# @pytest.mark.dragons_remote_data
# def test_tag_as_standard_real():
#     path = astrodata.testing.download_from_archive(
#         "xxx.fits", url="https://docs.google.com/uc?export=download"
#     )
#     ad = astrodata.open(path)
#     assert 'STANDARD' in ad.tags


# def test_ra_dec_from_text(astrofaker):
#     ad = astrofaker.create('SCORPIO', ['SPECT'],
#                            extra_keywords={'RA': '03:48:30.113',
#                                            'DEC': '+24:20:43.00',
#                                            'DATE-OBS': '2021-01-01T12:00:00.000'}
#                            )
#     assert ad.ra() == pytest.approx(57.12547083333333)
#     assert ad.dec() == pytest.approx(24.345277777777778)

#     # test bad RA/DEC values, just doing this for GMOS but it's testing the base
#     ad = astrofaker.create('SCORPIO', ['SPECT'],
#                            extra_keywords={'RA': 'Fail',
#                                            'DEC': 'Fail',
#                                            'DATE-OBS': '2021-01-01T12:00:00.000'}
#                            )
#     assert ad.ra() is None
#     assert ad.dec() is None


if __name__ == '__main__':

    pytest.main()
