#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_image_nearIR.py
# ------------------------------------------------------------------------------

from gempy.gemini import gemini_tools as gt

from geminidr.core import Image, Photometry
from .primitives_scorpio_nearIR import ScorpioNearIR
from . import parameters_scorpio_nearIR_image

from recipe_system.utils.decorators import parameter_override

import numpy as np
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioNearIRImage(ScorpioNearIR, Image, Photometry):
    """
    This class contains primitives that applies to all Scorpio near-IR
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'IMAGE', 'NIR'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioNearIRImage, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR_image)

