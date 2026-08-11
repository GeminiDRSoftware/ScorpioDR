#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_spect.py
# ------------------------------------------------------------------------------
import os
from importlib import import_module

from geminidr.core import Spect
from gempy.library import wavecal
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

    def _get_linelist(self, wave_model=None, *args, **kwargs):
        """
        Returns a list of wavelengths of the arc reference lines used by the
        primitive `determineWavelengthSolution()`, if the user parameter
        `linelist=None` (i.e., the default list is requested).

        Parameters
        ----------
        wave_model : astropy.modeling.models.Chebyshev1D instance
            model (with domain) defining the wavelength (range) required

        Returns
        -------
        gempy.library.wavecal.LineList object
            arc line wavelengths (and optional weights)
        """
        lookup_dir = os.path.dirname(import_module('.__init__',
                                                   self.inst_lookups).__file__)
        filename = os.path.join(lookup_dir, 'CuAr_GMOS.dat')
        return wavecal.LineList(filename)
