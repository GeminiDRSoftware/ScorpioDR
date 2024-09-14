#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_image_ccd.py
# ------------------------------------------------------------------------------
from contextlib import suppress

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
            datasec_kw = ad._keyword_for('data_section')
            if "CAL" in ad.tags:
                ad = super().trimOverscan([ad], suffix=suffix)[0]
                for ext in ad:
                    with suppress(AttributeError, KeyError):
                        del ext.hdr['OVRSECS*']
                        del ext.hdr['OVRSECP*']
            else:
                # KL: this is not taking into account the OBJMASK
                # (which core trimOverscan handles)  Matters for QAP if IQ measured
                # before overscan correction.
                for ext in ad:
                    nints = ext.data.shape[0]
                    data_list = []
                    mask_list = []
                    variance_list = []
                    for i in range(nints):
                        temp_ad = astrodata.create(ad.phu)
                        # KL: deepcopy of ext.nddata needed, otherwise datasec reset in
                        # super trim and following integrations use that new trimmed
                        # section instead of the original untrimmed.
                        temp_ad.append(deepcopy(ext.nddata[i]), header=deepcopy(ext.hdr))
                        os_trimmed = super().trimOverscan(adinputs=[temp_ad], suffix=suffix)
                        data_list.append(os_trimmed[0][0].data)
                        mask_list.append(os_trimmed[0][0].mask)
                        variance_list.append(os_trimmed[0][0].variance)
                        if i == 0:
                            datsec = os_trimmed[0].hdr.get(datasec_kw)[0]
                    ext.reset(np.array(data_list), mask=np.array(mask_list), variance=np.array(variance_list))

                    with suppress(AttributeError, KeyError):
                        del ext.hdr['OVRSECS*']
                        del ext.hdr['OVRSECP*']
                        # need to adjust DATSEC*  or deleted them.
                        del ext.hdr[datasec_kw+'*']

                    ext.hdr.set(datasec_kw, datsec, comment=self.keyword_comments.get(datasec_kw))
                    ext.hdr.set('TRIMSEC', datsec, comment=self.keyword_comments.get('TRIMSEC'))

                ad.phu.set('TRIMMED', 'yes', self.keyword_comments['TRIMMED'])
                gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
                ad.update_filename(suffix=suffix, strip=True)

        return adinputs
