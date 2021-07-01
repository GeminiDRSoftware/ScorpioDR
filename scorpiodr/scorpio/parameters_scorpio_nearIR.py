# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_nearIR.py file.

from gempy.library import config
from geminidr.gemini import parameters_gemini
#from geminidr.core import parameters_ccd  # import core pkgs as needed.

class standardizeStructureConfig(parameters_gemini.standardizeStructureConfig):
    def setDefaults(self):
        self.attach_mdf = False

class referencePixelCorrectConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_refpixelCorrected", optional=True)

class myNewPrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

