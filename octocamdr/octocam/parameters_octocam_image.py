# This parameter file contains the parameters related to the primitives
# defined in the primitives_octocam_image.py file.

from geminidr.core.parameters_image import ParametersImage
from geminidr.core.parameters_photometry import ParametersPhotometry
from .parameters_octocam import ParametersOCTOCAM

class ParametersOCTOCAMImage(ParametersOCTOCAM,
                             ParametersImage, ParametersPhotometry):
    pass
