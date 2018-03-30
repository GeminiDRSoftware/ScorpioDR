# This parameter file contains the parameters related to the primitives
# defined in the primitives_octocam_spect.py file.

from geminidr.core.parameters_spect import ParametersSpect
from .parameters_octocam import ParametersOCTOCAM

class ParametersOCTOCAMSpect(ParametersOCTOCAM, ParametersSpect):
    pass
