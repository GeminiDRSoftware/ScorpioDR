#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_image_ccd.py
# ------------------------------------------------------------------------------
from contextlib import suppress
from copy import deepcopy

import numpy as np

import astrodata
from geminidr.core import CCD
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio import Scorpio
from . import parameters_scorpio_ccd
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCD(Scorpio, CCD):
    """
    This class contains primitives that applies to all Scorpio optical
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'CCD'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioCCD, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_ccd)

    def biasCorrect(self, adinputs=None, suffix=None, bias=None, do_cal=None):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        if bias is None:
            bias_list = self.caldb.get_processed_bias(adinputs)
        else:
            bias_list = (bias, None)

        for ad, bias, origin in zip(*gt.make_lists(adinputs, *bias_list, force_ad=(1,))):
            ad = super().biasCorrect([ad], suffix=suffix, bias=bias, do_cal=do_cal)

        return adinputs

    def subtractOverscan(self, adinputs=None, **params):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        sfx = params["suffix"]

        for ad in adinputs:
            ad = super().subtractOverscan([ad], **params)[0]

        return adinputs

    def trimOverscan(self, adinputs=None, suffix=None):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            datasec_kw = ad._keyword_for('data_section')
            ad = super().trimOverscan([ad], suffix=suffix)[0]
            for ext in ad:
                with suppress(AttributeError, KeyError):
                    del ext.hdr['OVRSECS*']
                    del ext.hdr['OVRSECP*']
                    # The original DATSECn in untrimmed coords are superseded
                    # by DATSEC (the descriptor returns a single, consolidated
                    # section whether there is one keyword or four):
                    del ext.hdr[datasec_kw+'?']

        return adinputs
