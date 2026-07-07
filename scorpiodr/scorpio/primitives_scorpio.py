#
#                                                                       DRAGONS
#
#                                                         primitives_scorpio.py
# ------------------------------------------------------------------------------
import numpy as np

import astrodata
from geminidr.gemini.primitives_gemini import Gemini
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from . import parameters_scorpio
from .lookups import timestamp_keywords as scorpio_stamps
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

    def _initialize(self, adinputs, **kwargs):
        #super(Scorpio, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
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
            ad = super().darkCorrect([ad], suffix=suffix, dark=dark, do_cal=do_cal)

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
            ad = super().flatCorrect([ad], suffix=suffix, flat=flat, do_cal=do_cal)

        return adinputs

    def stackIntegrations(self, adinputs=None, **params):
        # This primitive should only be run on Scorpio science images.

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]
        sfx = params["suffix"]

        # Collect together & stack multiple integrations taken at the same
        # pointing (based on their data labels). Modified from the original
        # version that stacked cube planes when the format was 4-dimensonal.
        exp_groups = {}
        for ad in adinputs:
            try:
                expid, n_int = ad.data_label().rsplit('-', 1)
            except AttributeError:
                raise ValueError(f"Failed to read data label for {ad.filename}")

            if expid not in exp_groups:
                exp_groups[expid] = []

            exp_groups[expid].append(ad)

        adoutputs = []
        for expid, group in exp_groups.items():
            new_ad = self.stackFrames(group, **params)[0]
            new_ad.update_filename(suffix=sfx, strip=True)
            gt.mark_history(new_ad, primname=self.myself(), keyword=timestamp_key)
            adoutputs.append(new_ad)

        return adoutputs

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

        # The SCORPIO-specific code here, for handling higher dimensions, was
        # removed when the input data format changed from 4D to 2D in 2026.

        return adinputs

    def standardizeWCS(self, adinputs=None, **params):
        return adinputs

    def transferAttribute(self, adinputs=None, source=None, **params):
        """
        This primitive takes an attribute (e.g., "mask", or "OBJCAT") from
        the AD(s) in another ("source") stream and applies it to the ADs in
        this stream. There must be the same number of ADs in each stream, or
        only 1 in the source stream, or (where there are multiple
        sub-integrations per "exposure") exactly 1 per unique exposure in the
        source stream, with the primary stream grouped by exposure.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        source: str
            name of stream containing ADs whose attributes you want
        attribute: str
            attribute to transfer from ADs in other stream
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        if source not in self.streams.keys():
            log.info(f"Stream {source} does not exist so nothing to transfer")
            return adinputs

        source_length = len(self.streams[source])
        exp_groups = {}
        if not (source_length == 1 or source_length == len(adinputs)):
            fail = False
            for ad in adinputs:
                try:
                    expid, n_int = ad.data_label().rsplit('-', 1)
                except AttributeError:
                    log.warning(f"Failed to read data label for {ad.filename}")
                    fail = True
                    break

                if expid in exp_groups:
                    if expid != last_expid:
                        log.warning("Exposure groups out of order")
                        fail = True
                        break
                else:
                    exp_groups[expid] = []
                    last_expid = expid

                exp_groups[expid].append(ad)

            if fail or (source_length != len(exp_groups)):
                log.warning("Incompatible stream lengths/order: "
                            f"{len(adinputs)} and {source_length}")
                return adinputs

        if exp_groups:
            adoutputs = adinputs
            source_stream = self.streams[source]
            for expid, sourcead in zip(exp_groups, source_stream):
                self.streams['expgroup'] = exp_groups[expid]
                self.streams[source] = [sourcead]
                super().transferAttribute(stream='expgroup',
                                          source=source,
                                          **params)
            del self.streams['expgroup']
            self.streams[source] = source_stream
        else:
            adinputs = super().transferAttribute(adinputs, source, **params)

        return adinputs
