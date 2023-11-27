#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_image_ccd.py
# ------------------------------------------------------------------------------
from contextlib import suppress

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt

from geminidr.core import CCD
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_ccd
from copy import deepcopy

import astrodata
import numpy as np

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCD(Scorpio, CCD):
    """
    This class contains primitives that applies to all Scorpio optical
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'CCD'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioCCD, self).__init__(adinputs, **kwargs)
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
            if "CAL" in ad.tags:
                ad = super().biasCorrect([ad], suffix=suffix, bias=bias, do_cal=do_cal)
            else:
                for ext in ad:
                    nints = ext.data.shape[0]
                    data_list, mask_list, variance_list = [], [], []
                    for i in range(nints):
                        temp_ad = astrodata.create(ad.phu)
                        temp_ad.append(ext.nddata[i], header=ext.hdr)
                        bias_corrected = super().biasCorrect(adinputs=[temp_ad], suffix=suffix, bias=bias, do_cal=do_cal)
                        data_list.append(bias_corrected[0][0].data)
                        mask_list.append(bias_corrected[0][0].mask)
                        variance_list.append(bias_corrected[0][0].variance)
                    ext.reset(np.array(data_list), mask=np.array(mask_list), variance=np.array(variance_list))
                ad.phu.set('BIASIM', bias.filename, self.keyword_comments['BIASIM'])
                gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
                ad.update_filename(suffix=suffix, strip=True)
        return adinputs

    def subtractOverscan(self, adinputs=None, **params):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        sfx = params["suffix"]

        for ad in adinputs:
            if "CAL" in ad.tags:
                ad = super().subtractOverscan([ad], **params)[0]
            else:
                for ext in ad:
                    nints = ext.data.shape[0]
                    data_list = []
                    mask_list = []
                    variance_list = []
                    for i in range(nints):
                        temp_ad = astrodata.create(ad.phu)
                        temp_ad.append(ext.nddata[i], header=deepcopy(ext.hdr))
                        os_subtracted = super().subtractOverscan(adinputs=[temp_ad], **params)
                        data_list.append(os_subtracted[0][0].data)
                        mask_list.append(os_subtracted[0][0].mask)
                        variance_list.append(os_subtracted[0][0].variance)
                    ext.reset(np.array(data_list), mask=np.array(mask_list), variance=np.array(variance_list))
                gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
                ad.update_filename(suffix=sfx, strip=True)
        return adinputs

    def trimOverscan(self, adinputs=None, suffix=None):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if "CAL" in ad.tags:
                ad = super().trimOverscan([ad], suffix=suffix)[0]
                for ext in ad:
                    with suppress(AttributeError, KeyError):
                        del ext.hdr['OVRSECS*']
                        del ext.hdr['OVRSECP*']
            else:
                for ext in ad:
                    nints = ext.data.shape[0]
                    data_list = []
                    mask_list = []
                    variance_list = []
                    for i in range(nints):
                        temp_ad = astrodata.create(ad.phu)
                        temp_ad.append(ext.nddata[i], header=deepcopy(ext.hdr))
                        os_trimmed = super().trimOverscan(adinputs=[temp_ad], suffix=suffix)
                        data_list.append(os_trimmed[0][0].data)
                        mask_list.append(os_trimmed[0][0].mask)
                        variance_list.append(os_trimmed[0][0].variance)
                    ext.reset(np.array(data_list), mask=np.array(mask_list), variance=np.array(variance_list))

                    with suppress(AttributeError, KeyError):
                        del ext.hdr['OVRSECS*']
                        del ext.hdr['OVRSECP*']

                ad.phu.set('TRIMMED', 'yes', self.keyword_comments['TRIMMED'])
                gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
                ad.update_filename(suffix=suffix, strip=True)

        return adinputs

    def myNewPrimitive(self, adinputs=None, **params):
        """
        Description...

        Parameters
        ----------
        suffix: str
            suffix to be added to output files

        param2: blah
            blah, blah

        Returns
        -------

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get params out
        param2 = params['param2']

        # Initialize the list of output AstroData objects
        # It is also possible to modify adinputs in place.
        adoutputs = []

        for ad in adinputs:

            # Do whatever checks on the input are necessary, for example:
            # Check whether this primitive as been run already.
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has"
                            "already been processed by myNewPrimitive".
                            format(ad.filename))
                continue

            # -----------------------
            # DR algorithm goes here
            # -----------------------

            # Timestamp
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

            adoutputs.append(ad_out)

        return adoutputs