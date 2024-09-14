#
#                                                                       DRAGONS
#
#                                                         primitives_scorpio.py
# ------------------------------------------------------------------------------

import astrodata
import numpy as np
from geminidr.gemini.primitives_gemini import Gemini
from geminidr.core.primitives_visualize import Visualize
from gempy.adlibrary.manipulate_ad import remove_single_length_dimension
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

    def display(self, adinputs=None, **params):
        """
        Displays an image on the ds9 display, using multiple frames if
        there are multiple extensions or integrations. Saturated pixels can be
        displayed in red, and overlays can also be shown.

        SCORPIO raw data is 4D (the axes are [integrations, up-the-ramp, y, x]),
        so it needs to be reduced to 2D frames in their own extensions before
        passing it up to the DRAGONS `display` primitive for display. Each type
        of data (visible, IR) has its own complication to handle: the visible
        data needs to have the overscan subtracted for each quadrant, while the
        IR data has multiple up-the-ramp frames that need handling in order to
        display it.

        Parameters
        ----------
        ignore: bool
            setting to True turns off the display
        remove_bias: bool
            attempt to subtract bias before displaying?
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))

        # No-op if ignore=True
        if params["ignore"]:
            log.warning("display turned off per user request")
            return

        remove_bias = params['remove_bias']
        integ_num = params['integration']

        display_ads = []
        for ad in adinputs:
            # Create a new astrodata object to hold any modified data that
            # might be needed, depending on the processing level
            new_ad = astrodata.create(ad.phu)
            new_ad.filename = ad.filename

            for ext in ad:
                # Simulated data uses unsigned ints, which can underflow. Recast
                # to int for safety.
                if ext.data.dtype.kind == 'u':  # any unsigned integer
                    ext.data = ext.data.astype(dtype=int, casting='same_kind',
                                               copy=False)

                # Start by checking the shape of the data to see how raw/processed
                # it is, which determines what we need to do to it.
                if len(ext.data.shape) == 4:
                    # Data shape is (integrations, up-the-ramp, y, x)
                    if 'CCD' in ad.tags:
                        # Squeeze the data to remove the empty UTR axis
                        ext.operate(np.squeeze)

                    if 'NIR' in ad.tags:
                        # Create a new data structure to hold the data in,
                        # minus the up-the-ramp axis which goes away.
                        new_data = np.empty([ext.data.shape[0],
                                             ext.data.shape[2],
                                             ext.data.shape[3]],
                                            dtype=int)
                        for i in range(ext.data.shape[0]):
                            # Perform a quick-'n'-dirty 'last minus first'
                            # subtraction of the up-the-ramp axis
                            log.debug("Subtracting last up-the-ramp frame "
                                      f"from first in integration {i+1}")
                            new_data[i, :, :] = np.subtract(ext.data[i, -1, :, :],
                                                            ext.data[i, 0, :, :],
                                                            dtype=(int, int))
                        ext.data = new_data

                if len(ext.data.shape) == 3:
                    # Data shape is (integrations, y, x)
                    # Split up frames along the integrations axis into extensions
                    if integ_num and integ_num > ext.data.shape[0]:
                        raise ValueError(f"Only {ext.data.shape[0]} integrations "
                                         f"in file {ad.filename}, attempted to "
                                         f"display integration #{integ_num}.")
                    for i in range(ext.data.shape[0]):
                        if integ_num is None or integ_num == i+1:
                            log.fullinfo(f"Displaying integration {integ_num}")
                            new_ad.append(ext)
                            new_ad[-1].data = ext.data[i, :, :]


                if len(ext.data.shape) == 2:
                    # Data have been processed to 2D images per extension
                    if integ_num and integ_num > len(ad):
                        raise ValueError(f"Only {len(ad)} integrations "
                                         f"in file {ad.filename}, attempted to "
                                         f"display integration #{integ_num}.")
                    if integ_num is None or integ_num == ext.id:
                        new_ad.append(ext)

        display_ads.append(new_ad)

        for ad in display_ads:
            if 'CCD' in ad.tags and remove_bias:
                for ext in ad:
                    # Check if already overscan corrected, and continue if so
                    if (ad.phu.get('BIASIM') or ad.phu.get('DARKIM') or
                        ad.phu.get(self.timestamp_keys["subtractOverscan"])):
                        log.fullinfo("Bias level has already been removed from "
                                     "data; no approximate correction will be "
                                     "performed")
                        continue

                    # Handle each of the four quadrants making up the array;
                    # they each have strips along their two outside edges
                    # for the overscan correction.
                    overscans = ext.overscan_section()
                    for j in range(4):
                        sec_s = overscans['serial'][j]
                        sec_p = overscans['parallel'][j]
                        arrs = np.ravel(ext.data[sec_s.y1:sec_s.y2,
                                                 sec_s.x1:sec_s.x2])
                        arrp = np.ravel(ext.data[sec_p.y1:sec_p.y2,
                                                 sec_p.x1:sec_p.x2])
                        bias_level = int(np.median(np.concatenate([arrs, arrp])))
                        log.debug(f"Subtracting overscan of {int(bias_level)}"
                                  f" from quadrant {j+1} of ext {ext.id}")
                        # Collect limits from the overscan sections, then
                        # take the extrema to get the region to subtract from.
                        ys = {y for sec in [sec_s, sec_p] for y in sec[2:]}
                        xs = {x for sec in [sec_s, sec_p] for x in sec[:2]}
                        ext.data[min(ys):max(ys), min(xs):max(xs)] -= bias_level


        super().display(adinputs=display_ads, **params)

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

    def stackIntegrations(self, adinputs=None, **params):
        # This primitive should only be run on Scorpio science images.

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]

        # Each input AD should contain 3-D extension(s), the first axis being
        # the integration axis. We're going to loop over the ADs and split the
        # 3-D extension along the first axis, creating a new temp list of ADs.
        # Then call stackFrames in order to collapse these temp ADs into a single
        # 2-D extension belonging to the original AD.
        for ad in adinputs:
            for ext in ad:
                ndims = len(ext.data.shape)
                try:
                    assert ndims == 3
                except AssertionError:
                    if ndims > 3:
                        log.warning("Scorpio stackIntegrations - expected 3 dimensions and more than 3 found. No changes will be made to this extension.")
                    else:
                        log.warning("Scorpio stackIntegrations - expected 3 dimensions and less than 3 found. No changes will be made to this extension.")
                    continue
                int_list = []
                nints = ext.data.shape[0]
                for i in range(nints):
                    temp_ad = astrodata.create(ad.phu)
                    temp_ad.append(ext.nddata[i], header=ext.hdr)
                    int_list.append(temp_ad)
                flattened = self.stackFrames(int_list, **params)
                ext.reset(flattened[0][0].nddata)

            ad.update_filename(suffix=sfx, strip=True)

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
