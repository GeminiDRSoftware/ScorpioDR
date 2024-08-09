# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio_ccd_image.py file.

from gempy.library import config
from geminidr.gemini import parameters_qa

class measureIQConfig(parameters_qa.measureIQConfig):
    remove_bias = config.Field("Remove estimated bias level before displaying?",
                               bool, True)

