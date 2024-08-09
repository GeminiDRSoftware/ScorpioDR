#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_spect_ccd.py
# ------------------------------------------------------------------------------

from gempy.gemini import gemini_tools as gt

from geminidr.core import Spect
from .primitives_scorpio_ccd import ScorpioCCD
from . import parameters_scorpio_ccd_spect

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCDSpect(ScorpioCCD, Spect):
    """
    This class contains primitives that applies to all Scorpio optical
    spectroscopy data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'SPECT', 'CCD'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioCCDSpect, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_ccd_spect)

