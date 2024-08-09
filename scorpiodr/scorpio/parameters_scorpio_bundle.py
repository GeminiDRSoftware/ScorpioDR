# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_bundle.py file.

from gempy.library import config

"""
class somePrimitive(config.Config):
    suffix = config.Field("Filename suffix", str, "_suffix")
    param1 = config.Field("Param1", str, "default")
    param2 = config.Field("do param2?", bool, False)
"""

class separateChannelsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_unbundled", optional=True)