#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_image_ccd.py
# ------------------------------------------------------------------------------
from geminidr.core import Image, Photometry
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio_ccd import ScorpioCCD
from . import parameters_scorpio_ccd_image
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCDImage(ScorpioCCD, Image, Photometry):
    """
    This class contains primitives that applies to all Scorpio optical
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'IMAGE', 'CCD'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioCCDImage, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_ccd_image)
