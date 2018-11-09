# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio.py file.

from gempy.library import config
from geminidr.core import parameters_standardize

class validateDataConfig(parameters_standardize.validateDataConfig):
    def setDefaults(self):
        self.num_exts = [1] # add number of extension in raw data.

class somePrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

