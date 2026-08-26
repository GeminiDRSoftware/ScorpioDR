# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_spect.py file.

from geminidr.core import parameters_spect
from gempy.library import config

class calculateSensitivityConfig(parameters_spect.calculateSensitivityConfig):
    order = config.RangeField("Order of fitting function", int, 4, min=1)

