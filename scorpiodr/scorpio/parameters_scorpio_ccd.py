# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_ccd.py file.

from gempy.library import config
from geminidr.core import parameters_ccd
from geminidr.gemini import parameters_gemini

class standardizeStructureConfig(parameters_gemini.standardizeStructureConfig):
    def setDefaults(self):
        self.attach_mdf = False

class overscanCorrectConfig(parameters_ccd.overscanCorrectConfig):
    def setDefaults(self):
        self.function = "spline3"
        self.order = 0
        self.bias_type = "serial"


