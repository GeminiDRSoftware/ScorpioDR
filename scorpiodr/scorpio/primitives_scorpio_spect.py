#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_spect.py
# ------------------------------------------------------------------------------
from geminidr.core import Spect
from recipe_system.utils.decorators import parameter_override

from . import parameters_scorpio_spect
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioSpect(Spect):
    """
    This class contains primitives that applies to all Scorpio
    spectroscopy data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'SPECT'])

    def _initialize(self, adinputs, **kwargs):
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_spect)

    def standardizeWCS(self, adinputs=None, **params):
        """
        This primitive updates the WCS attribute of each NDAstroData extension
        in the input AstroData objects. For spectroscopic data, it means
        replacing an imaging WCS with an approximate spectroscopic WCS.

        Parameters
        ----------
        suffix: str/None
            suffix to be added to output files
        """
        super().standardizeWCS(adinputs, **params)
        for ad in adinputs:
            self._add_longslit_wcs(ad, pointing="center")
        return adinputs
