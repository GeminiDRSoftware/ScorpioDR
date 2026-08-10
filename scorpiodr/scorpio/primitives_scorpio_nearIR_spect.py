#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_spect_nearIR.py
# ------------------------------------------------------------------------------
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio_spect import ScorpioSpect
from .primitives_scorpio_nearIR import ScorpioNearIR
from . import parameters_scorpio_nearIR_spect
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioNearIRSpect(ScorpioSpect, ScorpioNearIR):
    """
    This class contains primitives that applies to all Scorpio near-IR
    spectroscopy data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'SPECT', 'NIR'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioNearIRSpect, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR_spect)

