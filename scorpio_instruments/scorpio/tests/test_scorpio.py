import datetime as dt
from types import NoneType

import numpy as np
import pytest

import astrodata
import astrodata.testing
from gemini_instruments.common import Section

import scorpio_instruments


SCORPIO_DESCRIPTORS_TYPES = [
    ('airmass', float),
    ('amp_read_area', [[str]]),
    ('ao_seeing', NoneType),
    ('array_name', [[str]]),
    ('array_section', [[Section]]),
    ('azimuth', float),
    ('calibration_key', str),
    ('camera', str),
    ('cass_rotator_pa', float),
    ('central_wavelength', float),
    ('coadds', int),
    ('data_label', str),
    ('data_section', [Section]),  # descriptor converts DATSEC1-4 to 1 section
    ('dec', float),
    ('decker', int | None),
    ('detector_name', str),
    ('detector_roi_setting', str),
    ('detector_rois_requested', NoneType),  # not needed (would be [Section])?
    ('detector_section', [Section]),
    ('detector_x_bin', int),
    ('detector_y_bin', int),
    ('detector_x_offset', float),
    ('detector_y_offset', float),
    ('disperser', str),
  # ('dispersion', [float]),
    ('dispersion_axis', [int]),
  # ('effective_wavelength', float),
    ('elevation', float),
    ('exposure_time', float),
  # ('filter_name', str),
  # ('focal_plane_mask', str),
    ('gain', [[float]]),
  # ('gain_setting', str),
    ('gcal_lamp', str | None),
    ('group_id', str),
    ('instrument', str),
    ('is_ao', bool),
    ('is_coadds_summed', bool),
  # ('local_time', dt.time),
  # ('lyot_stop', str),
    ('mdf_row_id', int | None),
  # ('nod_count', [int]),      # } currently undefined except for GMOS
  # ('nod_offsets', [float]),  # }
    ('nominal_atmospheric_extinction', float),
  # ('nominal_photometric_zeropoint', [float]),
    ('non_linear_level', int | float),  # should this be a len-1 list?
    ('object', str),
    ('observation_class', str),
  # ('observation_epoch', float),
    ('observation_id', str),
    ('observation_type', str),
  # ('overscan_section', [Section]),  # returns a dict; currently unsupported
    ('pixel_scale', float),
    ('program_id', str),
    ('pupil_mask', NoneType),
    ('qa_state', str),
    ('ra', float),
    ('raw_bg', int | None),  # }
    ('raw_cc', int | None),  # } can be 'UNKNOWN' if not set -> None
    ('raw_iq', int | None),  # }
    ('raw_wv', int | None),  # }
  # ('read_mode', str),
    ('read_noise', [[float]]),
  # ('read_speed_setting', str),
    ('requested_bg', int),
    ('requested_cc', int),
    ('requested_iq', int),
    ('requested_wv', int),
    ('saturation_level', int | float),  # should we enforce a list?
  # ('shuffle_pixels', int),  # currently undefined except for GMOS
  # ('slit', str),
    ('target_dec', float),
    ('target_ra', float),
    ('telescope', str),
    ('telescope_x_offset', float),
    ('telescope_y_offset', float),
    ('ut_date', dt.date),
    ('ut_datetime', dt.datetime),
    ('ut_time', dt.time),
    ('wavefront_sensor', str | None),
  # ('wavelength_band', str),
    ('wcs_dec', float),
    ('wcs_ra', float),
  # ('well_depth_setting', str),
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


def check_type(value, expected_type):
    """
    Here we deliberately recognize *only lists* of type(s) as a shorthand to
    indicate where descriptors return a list/tuple (or list of lists), since
    the types can themselves be containers (such as Section), or a tuple of
    possible types can be passed through to isinstance in the usual fashion,
    eg. "expected_type=[[(float, int, NoneType)]]" or "[[float | int | None]]".
    This test doesn't check type consistency within a list. Where a descriptor
    can return either a list or a single value, the corresponding tests would
    need separating out, to allow specifying different expected_types.
    """
    if isinstance(expected_type, list):
        try:
            iter(value)
            if isinstance(value, (str, bytes)):
                raise TypeError  # don't treat string as a container
        except TypeError:
            assert isinstance(value, list)  # fail with details
        else:
            expected_type = expected_type[0]
            for item in value:
                check_type(item, expected_type)
    else:
        assert isinstance(value, expected_type)


@pytest.mark.parametrize("descriptor,expected_type", SCORPIO_DESCRIPTORS_TYPES)
@pytest.mark.dragons_remote_data
def test_descriptor_matches_type(ad, descriptor, expected_type):
    value = getattr(ad, descriptor)()
    # print(descriptor, value)
    try:
        check_type(value, expected_type)
    except AssertionError as e:
        raise AssertionError(
            f"Assertion failed for file: {ad.filename}: {e}"
        ) from e


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
