# This parameter file contains the parameters related to the primitives
# defined in the primitives_octocam_spect_ccd.py file.

from geminidr.core.parameters_ccd import ParametersCCD
from .parameters_octocam_spect import ParametersOCTOCAMSpect

class ParametersOCTOCAMSpectCCD(ParametersOCTOCAMSpect, ParametersCCD):

    myNewPrimitive = {
        "suffix"        : "_newPrim",
        "param2"        : 5,
        "param3"        : None,
        "param4"        : True
    }

    # It is also possible to override defaults from inherited primitives.