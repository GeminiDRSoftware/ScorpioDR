#
#                                                                       DRAGONS
#
#                                                         primitives_scorpio.py
# ------------------------------------------------------------------------------

import astrodata
import numpy as np
from geminidr.gemini.primitives_gemini import Gemini
from gempy.gemini import gemini_tools as gt

from . import parameters_scorpio

from .lookups import timestamp_keywords as scorpio_stamps

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class Scorpio(Gemini):
    """
    This class inherits from the level above.  Any primitives specific
    to Scorpio and that applies to all channels can go here.  Primitives
    set specific to visible or nearIR, imaging or spectroscopy will go
    into subclasses inheriting Scorpio.
    """

    tagset = set()  # Cannot be assigned as a class

    def __init__(self, adinputs, **kwargs):
        super(Scorpio, self).__init__(adinputs, **kwargs)
        self.keyword_comments['DATSEC'] = self.keyword_comments['DATASEC']
        self._param_update(parameters_scorpio)
        # Add Scorpio-specific timestamp keywords
        self.timestamp_keys.update(scorpio_stamps.timestamp_keys)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be verified. This would be for the debundled data
        # and inherited and used by the ScorpioImage and ScorpioSpect class.
        return len(ad) in [1]

    def darkCorrect(self, adinputs=None, suffix=None, dark=None, do_cal=None):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        if dark is None:
            dark_list = self.caldb.get_processed_dark(adinputs)
        else:
            dark_list = (dark, None)

        for ad, dark, origin in zip(*gt.make_lists(adinputs, *dark_list, force_ad=(1,))):
            for ext in ad:
                nints = ext.data.shape[0]
                data_list, mask_list, variance_list = [], [], []
                for i in range(nints):
                    temp_ad = astrodata.create(ad.phu)
                    temp_ad.append(ext.nddata[i], header=ext.hdr)
                    dark_corrected = super().darkCorrect(adinputs=[temp_ad], suffix=suffix, dark=dark, do_cal=do_cal)
                    data_list.append(dark_corrected[0][0].data)
                    mask_list.append(dark_corrected[0][0].mask)
                    variance_list.append(dark_corrected[0][0].variance)
                ext.reset(np.array(data_list), 
                                   mask=np.array(mask_list), 
                                   variance=np.array(variance_list))
            ad.phu.set('DARKIM', dark.filename, self.keyword_comments['DARKIM'])
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=suffix, strip=True)
        return adinputs

    def flatCorrect(self, adinputs=None, suffix=None, flat=None, do_cal=None):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        if flat is None:
            flat_list = self.caldb.get_processed_flat(adinputs)
        else:
            flat_list = (flat, None)

        for ad, flat, origin in zip(*gt.make_lists(adinputs, *flat_list, force_ad=(1,))):
            for ext in ad:
                nints = ext.data.shape[0]
                data_list, mask_list, variance_list = [], [], []
                for i in range(nints):
                    temp_ad = astrodata.create(ad.phu)
                    temp_ad.append(ext.nddata[i], header=ext.hdr)
                    flat_corrected = super().flatCorrect(adinputs=[temp_ad], suffix=suffix, flat=flat, do_cal=do_cal)
                    data_list.append(flat_corrected[0][0].data)
                    mask_list.append(flat_corrected[0][0].mask)
                    variance_list.append(flat_corrected[0][0].variance)
                ext.reset(np.array(data_list), 
                                   mask=np.array(mask_list), 
                                   variance=np.array(variance_list))
            ad.phu.set("FLATIM", flat.filename, self.keyword_comments["FLATIM"])
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=suffix, strip=True)
        return adinputs

    def standardizeInstrumentHeaders(self, adinputs=None, suffix=None):
        """
        This primitive is used to make the changes and additions to the keywords in the headers of SCORPIO data, specifically.

        2023-01-30 - Landon Gelman
        The contents of this function are currently intended just for testing SCORPIO's implementation with DRAGONS' core primitives.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """

        # Instantiate the log
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already  been processed by "
                            "standardizeInstrumentHeaders".format(ad.filename))
                continue

            # Standardize the headers of the input AstroData object. Update the 
            # keywords in the headers that are specific to SCORPIO.
            log.status("Updating keywords that are specific to SCORPIO")

            for desc in ('saturation_level', 'non_linear_level'):
                kw = ad._keyword_for(desc)
                comment = self.keyword_comments[kw]
                dv = getattr(ad, desc)()
                if isinstance(dv, list):
                    for ext, value in zip(ad, dv):
                        ext.hdr.set(kw, value, comment)
                else:
                    ad.hdr.set(kw, dv, comment)

            # Timestamp and update filename
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=suffix, strip=True)
        return adinputs

    def standardizeStructure(self, adinputs=None, **params):
        """
        This primitive is used to standardize the structure of SCORPIO data,
        specifically.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        attach_mdf: bool
            attach an MDF to the AD objects?
        mdf: str
            full path of the MDF to attach
        """
        adinputs = super().standardizeStructure(adinputs, **params)

        adoutputs = []
        for ad in adinputs:
            if "CCD" in ad.tags:     # VIS
                new_ad = astrodata.create(ad.phu)
                for ext in ad:
                    nints = ext.data.shape[0]
                    if "CAL" in ad.tags:
                        # Remove the group axis and integrations axis
                        # For CALs, there's only 1 integration per exposure and we need a stackable 2-D image
                        new_ad.append(ext.data[0, 0], header=ext.hdr)
                    else:
                        # Skip the dark frames preceding light frames in science frames in continuous window imaging.
                        if "IMAGE" in ad.tags and ad.phu["WMODE"].upper() == "CONTINUOUS":
                            new_ad.append([ext.data[nint, 0] for nint in range(ext.hdr["WNFDARK"], nints)], header=ext.hdr)
                        else:   # full frame imaging, burst window imaging, or spectroscopy
                            new_ad.append([ext.data[nint, 0] for nint in range(nints)], header=ext.hdr)
                new_ad.filename = ad.filename
                adoutputs.append(new_ad)
            else:   # NIR
                adoutputs.append(ad)
        return adoutputs

    def standardizeWCS(self, adinputs=None, **params):
        return adinputs