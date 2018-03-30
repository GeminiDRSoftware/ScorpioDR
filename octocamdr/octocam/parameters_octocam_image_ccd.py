# This parameter file contains the parameters related to the primitives
# defined in the primitives_octocam_image_ccd.py file.

from geminidr.core.parameters_ccd import ParametersCCD
from .parameters_octocam_image import ParametersOCTOCAMImage

class ParametersOCTOCAMImageCCD(ParametersOCTOCAMImage, ParametersCCD):

    myNewPrimitive = {
        "suffix"        : "_newPrim",
        "param2"        : 5,
        "param3"        : None,
        "param4"        : True
    }

    # It is also possible to override defaults from inherited primitives.