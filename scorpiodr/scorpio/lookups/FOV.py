# This code is looked up by gempy as part of the configuration for the
# appropriate instrument and evaled by the infrastructure. It has initially
# been written to support gemini_tools.ExposureGroup and is semi-optimized.
#
# All calculations are now done in PIXELS

def pointing_in_field(ad, refpos, frac_FOV=1.0, frac_slit=1.0):

    """
    See gemini_tools.pointing_in_field() for the API. This is an
    instrument-specific back end that you shouldn't be calling directly.

    No inputs are validated at this level; that's the responsibility of the
    calling function, for reasons of efficiency.

    :type pos: AstroData instance
    :type refpos: tuple of floats
    :type frac_FOV: float
    :type frac_slit: float
    """

    # Extract pointing info from the AstroData instance.
    position = (ad.detector_x_offset(), ad.detector_y_offset())

    #dist = frac_FOV * 1024 # pixels

    return all(abs(x-r) < frac_FOV * ad[0].shape[0] 
               for x, r in zip(position, refpos))