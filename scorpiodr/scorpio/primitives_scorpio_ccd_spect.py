#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_spect_ccd.py
# ------------------------------------------------------------------------------
from geminidr.core import Spect
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio_ccd import ScorpioCCD
from . import parameters_scorpio_ccd_spect
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCDSpect(ScorpioCCD, Spect):
    """
    This class contains primitives that applies to all Scorpio optical
    spectroscopy data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'SPECT', 'CCD'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioCCDSpect, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_ccd_spect)

