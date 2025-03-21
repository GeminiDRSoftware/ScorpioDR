#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_image_nearIR.py
# ------------------------------------------------------------------------------
import numpy as np

from geminidr.core import Image, Photometry
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio_nearIR import ScorpioNearIR
from . import parameters_scorpio_nearIR_image
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioNearIRImage(ScorpioNearIR, Image, Photometry):
    """
    This class contains primitives that applies to all Scorpio near-IR
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'IMAGE', 'NIR'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioNearIRImage, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR_image)

