# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_nearIR.py file.

from gempy.library import config
from geminidr.gemini import parameters_gemini

class calculateSignalByRegressionConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_slopeDetermined", optional=True)

class flagCosmicRaysFromNDRsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_CRsFlagged", optional=True)

class standardizeStructureConfig(parameters_gemini.standardizeStructureConfig):
    def setDefaults(self):
        self.attach_mdf = False

class subtractReferencePixelsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_refpixelsSubtracted", optional=True)
    do_vertical_correction = config.Field("Perform vertical reference pixel corrections?", bool, True)

class trimReferencePixelsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_refpixelsTrimmed", optional=True)

class referencePixelsCorrectConfig(subtractReferencePixelsConfig, trimReferencePixelsConfig):
    suffix = config.Field("Filename suffix", str, "_refpixelsCorrected", optional=True)



