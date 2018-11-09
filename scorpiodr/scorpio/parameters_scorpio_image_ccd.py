# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_image_ccd.py file.

from gempy.library import config
#from geminidr.core import parameters_ccd  # import core pkgs as needed.

class myNewPrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)

