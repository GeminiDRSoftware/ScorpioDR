# This parameter file contains the parameters related to the primitives
# defined in the primitives_octocam_image.py file.

from gempy.library import config
from geminidr.core import parameters_standardize

class standardizeStructureConfig(parameters_standardize.standardizeStructureConfig):
    def setDefaults(self):
        self.attach_mdf = False

class myNewPrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

