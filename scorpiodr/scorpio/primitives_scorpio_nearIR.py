#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_spect_nearIR.py
# ------------------------------------------------------------------------------
import os
import sys
from copy import deepcopy
from datetime import timedelta

import numpy as np

from astropy.table import Table

import astrodata
from astrodata.nddata import NDAstroData, ADVarianceUncertainty
from geminidr.core import NearIR
from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import capture_provenance
from recipe_system.utils.decorators import parameter_override

from . import fitramp

from .primitives_scorpio import Scorpio
from . import parameters_scorpio_nearIR
# ------------------------------------------------------------------------------

@parameter_override
@capture_provenance
class ScorpioNearIR(Scorpio, NearIR):
    """
    This class contains primitives that applies to all Scorpio near-IR data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'NIR'])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioNearIR, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR)

    def associateSky(self, adinputs=None, **params):

        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        sky = params.pop('sky', [])
        if sky:
            # Produce a list of AD objects from the sky frame/list
            ad_skies = sky if isinstance(sky, list) else [sky]
            ad_skies = [ad if isinstance(ad, astrodata.AstroData) else
                           astrodata.open(ad) for ad in ad_skies]
        else:  # get from sky stream (put there by separateSky)
            ad_skies = self.streams.get('sky', [])

        # Timestamp and update filenames. Do now so filenames agree at end
        for ad in set(adinputs + ad_skies):
            ad.update_filename(suffix=params['suffix'], strip=True)
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

        # Split integrations stacked in the 3rd array dimension into separate
        # 2D images, as expected by the core primitive:
        split_obj, obj_map = self._split_integrations_to_ad(adinputs)
        split_sky, sky_map = self._split_integrations_to_ad(ad_skies)

        # Now run the core primitive on the split integrations, as usual:
        split_obj = super().associateSky(split_obj, sky=split_sky, **params)

        # Restore the original sky stream, since the core primitive overrides
        # it from the split_sky argument above:
        self.streams['sky'] = ad_skies

        # Combine the SKYTABLE entries for individual integrations into a set
        # for each input ad, with the integration numbers in separate columns:
        sky_data = {ad.filename : [] for ad in adinputs}
        for ad in split_obj:
            obj_ad, obj_int = obj_map[ad.filename]
            for sky_fn in ad.SKYTABLE['SKYNAME']:
                sky_ad, sky_int = sky_map[sky_fn]
                sky_data[obj_ad.filename].append(
                    (np.int16(obj_int), sky_ad.filename, np.int16(sky_int))
                )

        for ad in adinputs:
            ad.SKYTABLE = Table(
                names=('OBJINT', 'SKYNAME', 'SKYINT'),
                rows=sky_data[ad.filename]
            )

        return adinputs

    def calculateSignalByRegression(self, adinputs=None, **params):
        """
        For each pixel, iteratively fit a slope in electrons/s to a ramp of
        non-destructive read group averages, allowing for jumps due to cosmic
        rays, and convert to a number of electrons during the integration.

        See Brandt (2024), PASP 136, 045004 & 045005.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]

        # Consider moving this into a descriptor later:
        def _get_read_times(ext):
            n_groups = ext.shape[-3] if len(ext.shape) > 2 else 1
            if getattr(ext.hdr, 'UTRGROUP', n_groups) != n_groups:
                raise ValueError("Number of read groups UTRGROUP doesn't match "
                                 "array dimensions")
            read_interval = 0.38  # except in high-speed imaging (without UTR?)
            n_reads = ext.hdr['UTRFRAME']
            n_total = n_reads + ext.hdr['UTRSKIP']

            read_times = []
            elapsed = 0.
            for ngroup in range(n_groups):
                read_times.append([])
                for nread in range(n_total):
                    elapsed += read_interval
                    if nread < n_reads:
                        read_times[-1].append(elapsed)

            return read_times

        adoutputs = []
        for ad in adinputs:

            new_ad = astrodata.create(ad.phu)

            for ext in ad:

                # Don't convert ADU automatically here, to avoid replicating
                # the fiddling with headers & history in ADUToElectrons:
                if ext.is_in_adu():
                    raise ValueError(f'Please convert {ad.filename} to units '
                                     f'of electrons')

                n_ints, n_groups, n_rows, n_cols = ext.data.shape
                read_times = _get_read_times(ext)  # per integration
                n_reads = len(read_times[0])

                # Construct covariance matrix components needed by fit_ramps,
                # based on the sampling in time:
                covar_obj = fitramp.Covar(read_times)
                tgroup = covar_obj.delta_t[:, np.newaxis, np.newaxis]

                # Get array of read noise by pix (amp) and convert from 4D-2D:
                rdnoise = gt.array_from_descriptor_value(ext, "read_noise")
                rdnoise = rdnoise[tuple(0 for dim in rdnoise.shape[:-2])]

                # TO DO: CONFIRM WHETHER READ NOISE VALUES CORRESPOND TO THE
                # GROUP AVERAGE (as per next line) OR AN INDIVIDUAL READ.
                rdnoise *= np.sqrt(n_reads)

                # Create a new extension to hold the output data:
                outshape = (n_ints, n_rows, n_cols)
                new_ad.append(
                    NDAstroData(
                        data = np.empty(outshape, dtype=ext.data.dtype),
                        mask = (None if ext.mask is None else
                                np.zeros(outshape, dtype=ext.mask.dtype)),
                        variance = np.zeros(outshape, dtype=ext.variance.dtype)
                    ),
                    header=ext.hdr
                )
                new_ext = new_ad[-1]

                # Remove the NDR group axis from the WCS in the output:
                new_ext.wcs = ext.nddata.window[:, 0].wcs
                astrodata.wcs.remove_unused_world_axis(new_ext)

                # Keep stats on goodness of fits for the log:
                chisq = np.empty((n_rows, n_cols), dtype=np.float32)

                for n_int, (integ, new_integ) in enumerate(
                        zip(ext.nddata, new_ext.nddata)
                ):

                    log.stdinfo(f'Integration {n_int}')

                    diff_rates = ((integ.data[1:] - integ.data[:-1]) / tgroup)

                    # Track which diffs are used in the final fit, masking any
                    # beforehand that come from saturated values (1 = good),
                    # unless all of them are saturated for a given pixel:
                    if integ.mask is None:
                        diff_mask = np.ones_like(diff_rates, dtype=np.uint8)
                    else:
                        diff_dq = integ.mask[1:] | integ.mask[:-1]
                        diff_mask = diff_dq & DQ.saturated
                        diff_mask = (
                            np.logical_not(diff_mask) |
                            np.all(diff_mask, axis=0)[np.newaxis, ...]
                        ).astype(np.uint8)

                    for nrow in range(n_rows):  # TO DO: transpose arrays?

                        # This gets modified by reference below:
                        good_diffs = diff_mask[:, nrow]

                        # Initial ramp fitting, allowing for jumps due to
                        # cosmic rays (by comparing the goodness of fit when
                        # excluding diffs between consecutive read groups):
                        good_diffs, fitted_rate = fitramp.mask_jumps(
                            diff_rates[:, nrow], covar_obj, rdnoise[nrow],
                            diffs2use=good_diffs
                        )
                        # Next iter requires non-negative count rate estimates:
                        fitted_rate = np.clip(fitted_rate, a_min=0, a_max=None)

                        # Final iteration (with CR jumps already masked), which
                        # produces noise/X^2 and corrects for statistical bias:
                        result = fitramp.fit_ramps(
                            diff_rates[:, nrow], covar_obj, rdnoise[nrow],
                            diffs2use=good_diffs,
                            countrateguess=fitted_rate
                        )

                        # TO DO: Convert back to e- with exp_time descriptor:
                        new_integ.data[nrow] = result.countrate
                        new_integ.variance[nrow] = result.uncert**2

                        # Combine DQ from all the diffs actually used in each
                        # ramp fit, after excluding saturated reads & CR jumps:
                        if integ.mask is not None:
                            new_integ.mask[nrow] = np.bitwise_or.reduce(
                                diff_dq[:, nrow] * good_diffs,
                                axis=0
                            )

                        chisq[nrow] = result.chisq

                    # Report chi^2 divided by by the DOF of each fit, ignoring
                    # any pixels that don't have multiple good diffs:
                    dof = diff_mask.sum(axis=0) - 1
                    good_diffs = dof > 0  # (note: briefly re-using variable)
                    rchisq = chisq[good_diffs] / dof[good_diffs]
                    log.stdinfo(f'  red. chi^2: {rchisq.mean():.2f} +/-'
                                f' {rchisq.std():.2f} (range'
                                f' {rchisq.min():.2f}--{rchisq.max():.1f})')
                    del rchisq

                # For calibrations only, check that the integration axis is
                # redundant and remove it (as well as the read group axis),
                # leaving a 2D array. The astrodata function below can only
                # remove one axis at a time.
                # Eventually this should be moved elsewhere in ScorpioDR?
                if "CAL" in ad.tags:
                    try:
                        assert n_ints == 1
                    except AssertionError:
                        msg = "CAL integration axis dimension exceeds 1"
                        log.error(msg)
                        raise ValueError(msg)
                    new_ext.reset(new_ext.nddata[0])
                    astrodata.wcs.remove_unused_world_axis(new_ext)

            # Update the filename.
            ad.update_filename(suffix=sfx, strip=True)
            new_ad.filename = ad.filename
            adoutputs.append(new_ad)

        return adoutputs

    def referencePixelsCorrect(self, adinputs=None, **params):
        adinputs = self.subtractReferencePixels(adinputs, 
                    **self._inherit_params(params, "subtractReferencePixels"))
        adinputs = self.trimReferencePixels(adinputs, suffix=params["suffix"])
        return adinputs

    def skyCorrect(self, adinputs=None, **params):
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        # timestamp_key = self.timestamp_keys[self.myself()]

        suffix = params["suffix"]

        ad_skies = self.streams["sky"] if "sky" in self.streams else adinputs

        # Split integrations stacked in the 3rd array dimension into separate
        # 2D images, as expected by the core primitive:
        split_obj, obj_map = self._split_integrations_to_ad(adinputs)
        split_sky, sky_map = self._split_integrations_to_ad(ad_skies)

        self.streams["sky"] = split_sky  # temp hack till sky param fixed

        split_obj = super().skyCorrect(split_obj, **params)

        self.streams["sky"] = ad_skies   # still needed anyway?

        # Update the mapping keys with the new file suffix, after running the
        # core skyCorrect (don't just update names before splitting, because
        # skyCorrect expects specific suffixes):
        del sky_map  # just so it's not inconsistent
        obj_map = {
            self._modify_filename_components(
                split_fn, n_int=None, suffix=suffix
            ) : entry for split_fn, entry in obj_map.items()
        }

        # Although the split ad instances are views into the original adinputs,
        # the sky subtraction can make copies and appears not to change the
        # parent objects (which could make a brittle assumption anyway), so
        # copy the results (& some meta-data updates) into adinputs explicitly:
        for ad in split_obj:
            orig_ad, n_int = obj_map[ad.filename]
            if n_int == 0:
                orig_ad.update_filename(suffix=suffix, strip=True)
                orig_ad.phu.update({
                    k : (ad.phu[k], ad.phu.comments[k])
                    for k in set(ad.phu)-set(orig_ad.phu)
                })  # replaces "mark_history()"
            if len(ad) != len(orig_ad):
                raise ValueError(  # shouldn't happen, but just for clarity
                    "number of integrations must not vary between extensions"
                )
            for split_ext, orig_ext in zip(ad, orig_ad):
                idx = n_int if len(orig_ext.shape) > 2 else slice(None)
                # NDData apparently doesn't support assignment to a section!
                orig_ext.data[idx] = split_ext.data
                if orig_ext.variance is not None:
                    orig_ext.variance[idx] = split_ext.variance
                if orig_ext.mask is not None:
                    orig_ext.mask[idx] = split_ext.mask

        return adinputs

    def subtractReferencePixels(self, adinputs=None, **params):
        """
        Correct a SCORPIO NIR image's noise by using its reference pixels. 

        This algorithm takes some extra steps that wouldn't be necessary if the user
        is expecting to trim the reference pixels (and the unilluminated amplifiers in
        full frame imaging) immediately after this primitive is run. But in the event
        that a user wants or needs to inspect the image outputted from this primitive,
        the extra steps ensure that all pixels are corrected.

        The methodology used in this primitive is based on Robberto 2014,
        "On the Reference Pixel Correction of NIRCam Detectors" and his IDL
        code, corrector.pro.

        Robberto 2014:
            https://www.stsci.edu/files/live/sites/www/files/home/jwst/documentation/technical-documents/_documents/JWST-STScI-003852.pdf
        corrector.pro IDL code: 
            https://www.stsci.edu/~robberto/Main/Software/IDL4pipeline/

        Note, SCORPIO's reference pixels are laid out very specifically. Even
        the "full" image size is smaller than the size of the detector.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]
        do_vertical_correction = params["do_vertical_correction"]

        # Construct a 2D "time" array that matches the size and shape of one amplifier. 
        # Multiply it by 10E-6 (seconds) to approximate the pixel readout rate.
        # Slice out the corresponding horizontal and vertical reference pixels.
        # Calculate the row averages for the corresponding reference pixels.
        # Assemble the single column of time pixels as they would appear in the image, i.e. keepdims.
        # So we have an array of 2048 rows by 1 column.
        # Finally, get a 1D array of all pixels corresponding to the horizontal reference pixels for linear regression later on.
        amp_rows = 2048
        amp_cols = 64
        amp_size = amp_rows * amp_cols
        time = np.flip(np.linspace(1., amp_size, amp_size, endpoint=True).reshape((amp_rows, amp_cols)) * 10E-6, 0)
        time_top = time[:4, :]
        time_bot = time[-4:, :]
        time_side = time[4:-4, :4]
        time_top_row_avg = np.mean(time_top, axis=1, keepdims=True)
        time_bot_row_avg = np.mean(time_bot, axis=1, keepdims=True)
        time_side_row_avg = np.mean(time_side, axis=1, keepdims=True)
        time_row_avg = np.vstack([time_top_row_avg, time_side_row_avg, time_bot_row_avg])
        time_horizontal = np.vstack([time_top, time_bot]).flatten()    # used later for linear regression

        for ad in adinputs:
            obsmode = ad.phu["OBSMODE"].upper()    # imaging vs spectroscopy
            if obsmode == "IMAGE":
                imtype = ad.phu["IMTYPE"].upper()    # full frame vs window

            # Get the list of Section tuples for all reference pixels sections.
            refpix_sec = ad.refpix_section()

            for e, ext in enumerate(ad):
                # Make an empty_like array which will be used to replace the original extension data
                new_ext_data = np.empty_like(ext.data, dtype=np.float32)

                # Make a copy of the data so we don't modify the original yet.
                if ext.is_in_adu():
                    data = deepcopy(ext.data.astype(np.float32))    # 4D numpy array
                else:
                    data = deepcopy(ext.data)

                # Convert any saturated or bad pixels to NaNs so they are not used in any calculations.
                data[np.where(np.bitwise_and(ext.mask, DQ.saturated))] = np.nan
                data[np.where(np.bitwise_and(ext.mask, DQ.bad_pixel))] = np.nan

                # Get the indices of the reference pixel sections for this extension from the master list above.
                rpix_top_sec = refpix_sec["top"][e]     # top horiztonal ref pixels    # list of Section tuples
                rpix_bot_sec = refpix_sec["bottom"][e]     # bottom horizontal ref pixels

                # if full frame imaging or spectroscopy and vertical ref pix corrections on
                if ((obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT") and do_vertical_correction:
                    rpix_side_sec = refpix_sec["side"][e]    # vertical ref pixels

                for intn in range(len(ext.data)):
                    # Extract the reference pixel section from the data. These are lists of 3D arrays with the shape (ngroups, nrows, ncols).
                    # The arrays in rpt and rpb correspond to amplifiers. The arrays in rps do not yet correspond to amplifiers.
                    # Calculate and subtract the average for each pixel along the group axis.
                    rpt = [data[intn, :, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_top_sec]    # list of 3D numpy arrays
                    rpb = [data[intn, :, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_bot_sec]

                    ngroups = len(data[intn])
                    namps = len(rpt)

                    rpt_no_ga, rpb_no_ga = [], []    # group averaged horizontal reference pixels
                    for amp in range(namps):
                        meant = np.mean(rpt[amp], axis=0)
                        meanb = np.mean(rpb[amp], axis=0)
                        rpt_no_ga.append(rpt[amp]-meant)
                        rpb_no_ga.append(rpb[amp]-meanb)

                    # if full frame imaging or spectroscopy and vertical ref pix corrections on
                    if ((obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT") and do_vertical_correction:
                        rps_raw = [data[intn, :, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_side_sec]

                        # Restructure the vertical reference pixel sections to correspond to amplifiers.
                        rpix_left, rpix_right = [], []
                        for sec in range (len(rps_raw)):
                            if sec % 2 == 0:
                                rpix_right.append(rps_raw[sec])
                            else:
                                rpix_left.append(rps_raw[sec])
                        rps = [np.hstack(rpix_left), np.hstack(rpix_right)]

                        rps_no_ga = []    # group averaged vertical reference pixels
                        meanl = np.mean(rps[0], axis=0)
                        meanr = np.mean(rps[1], axis=0)
                        rps_no_ga.append(rps[0]-meanl)
                        rps_no_ga.append(rps[1]-meanr)

                    corrector_cube = np.empty_like(data[intn])
                    # Loop over groups
                    for grp in range(ngroups):
                        # Calculate the linear regression coefficients.
                        coeffs = []
                        for amp in range(namps):
                            rpix_horizontal = np.vstack([rpt_no_ga[amp][grp], rpb_no_ga[amp][grp]]).flatten()
                            coeffs.append(np.polyfit(time_horizontal, rpix_horizontal, deg=1))

                        # if full frame imaging or spectroscopy and vertical ref pix corrections on
                        if ((obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT") and do_vertical_correction:
                            # Calculate the row averges for the top, side, and bottom reference pixels
                            # in the first array.
                            row_avg_top_left = np.mean(rpt_no_ga[0][grp], axis=1, keepdims=True)
                            row_avg_bot_left = np.mean(rpb_no_ga[0][grp], axis=1, keepdims=True)
                            row_avg_side_left = np.mean(rps_no_ga[0][grp], axis=1, keepdims=True)
                            row_avg_left = np.vstack([row_avg_top_left, row_avg_side_left, row_avg_bot_left])

                            # Repeat for the last amplifier
                            row_avg_top_right = np.mean(rpt_no_ga[1][grp], axis=1, keepdims=True)
                            row_avg_bot_right = np.mean(rpb_no_ga[1][grp], axis=1, keepdims=True)
                            row_avg_side_right = np.mean(rps_no_ga[1][grp], axis=1, keepdims=True)
                            row_avg_right = np.vstack([row_avg_top_right, row_avg_side_right, row_avg_bot_right])

                            # Construct a "high frequency" function for the first and last amplifiers
                            # and average them row-wise.
                            highfreq_left = row_avg_left - coeffs[0][1] - (coeffs[0][0] * time_row_avg)
                            highfreq_right = row_avg_right - coeffs[-1][1] - (coeffs[-1][0] * time_row_avg)
                            highfreq = (highfreq_left + highfreq_right) / 2

                            # Now run the FFT smoothing function
                            delt = 10E-6 * amp_cols
                            smoothed_data = self._smoothFFT(highfreq, delt).reshape(highfreq.shape)

                        # Construct the 2D corrector arrays by amplifier
                        corrector_by_amp = []
                        for amp in range(namps):
                            # if full frame iamging or spectroscopy and vertical ref pix corrections on
                            if ((obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT") and do_vertical_correction:
                                corr = coeffs[amp][1] + (coeffs[amp][0] * time_row_avg) + smoothed_data
                            else:
                                corr = coeffs[amp][1] + (coeffs[amp][0] * time_row_avg)
                            corr = np.repeat(corr, amp_cols, axis=1)
                            corrector_by_amp.append(corr)
                        corrector_group = np.hstack(corrector_by_amp)    # 2048 x namps

                        # Extract the portions of the corrector array that are covered by the observation mode.
                        # if full frame imaging or spectroscopy
                        if (obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT":
                            # See the section descriptor document for a picture of how this looks
                            top_left = corrector_group[4:512, :4]          # vertical ref pixels
                            bot_left = corrector_group[1536:2044, :4]      # vertical ref pixels
                            top_right = corrector_group[4:512, -4:]        # vertical ref pixels
                            bot_right = corrector_group[1536:2044, -4:]    # vertical ref pixels

                            top_zero = np.zeros((4, 4)) # zero padding top vertical ref pixels
                            mid_zero = np.zeros((8, 4)) # zero padding middle vertical ref pixels
                            bot_zero = np.zeros((4, 4)) # zero padding bottome vertical ref pixels

                            left_ref = np.concatenate([top_zero, top_left, mid_zero, bot_left, bot_zero], axis=0)
                            right_ref = np.concatenate([top_zero, top_right, mid_zero, bot_right, bot_zero], axis=0)

                            top_ref = corrector_group[:4, :]     # horizontal ref pixels
                            bot_ref = corrector_group[-4:, :]    # horizontal ref pixels

                            sci = corrector_group[512:1536, :]    # active pixels and middle vertical ref pixels

                            corrector_group_cropped = np.concatenate([top_ref, sci, bot_ref], axis=0)
                            corrector_group_cropped = np.concatenate([left_ref, corrector_group_cropped, right_ref], axis=1)

                        # window imaging mode
                        else:
                            top_ref = corrector_group[:4, :]
                            bot_ref = corrector_group[-4:, :]
                            sci = corrector_group[970:1078, :]
                            corrector_group_cropped = np.concatenate([top_ref, sci, bot_ref], axis=0)

                        try:
                            assert ext.data[intn][grp].shape == corrector_group_cropped.shape
                        except AssertionError:
                            log.error("Data shape does not match corrector shape")
                            raise ValueError("Data shape does not match corrector shape")

                        corrector_cube[grp] = corrector_group_cropped

                    # Take the original data array and subtract the corrector cube.
                    # Then extract all horizontal reference pixels and get their average --> "superaverage"
                    # Then loop over the groups, and the amps, and get the average per amplifier per group.
                    # Calculate the "delta" subtracting the amp average from the superaverage.
                    data_corrected = data[intn] - corrector_cube
                    rpt_corrected = [data_corrected[:, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_top_sec]    # list of 3D numpy arrays
                    rpb_corrected = [data_corrected[:, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_bot_sec]
                    superaverage = np.mean([rpt_corrected, rpb_corrected])    # scalar
                    for grp in range(ngroups):
                        for amp in range(namps):
                            amp_avg = np.mean([rpt_corrected[amp][grp], rpb_corrected[amp][grp]])    # scalar
                            amp_delta = superaverage - amp_avg    # scalar

                            data_corrected[grp, :, rpix_top_sec[amp].x1:rpix_top_sec[amp].x2] += amp_delta

                            # if full frame imaging or spectroscopy
                            if (obsmode == "IMAGE" and imtype == "FULL") or obsmode == "SPECT":
                                if amp == 0:    # first amp
                                    data_corrected[grp, 4:512, :4] += amp_delta
                                    data_corrected[grp, 1536:2044, :4] += amp_delta
                                elif amp == namps-1:
                                    data_corrected[grp, 4:512, -4:] += amp_delta
                                    data_corrected[grp, 1536:2044, -4:] += amp_delta

                    new_ext_data[intn] = data_corrected

                ext.reset(new_ext_data)

            # Update filename
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs

    def trimReferencePixels(self, adinputs=None, **params):
        """
        This primitive is used to remove both the reference pixels and the
        non-illuminated active data pixels that are not included in the
        data section. 

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]

        for ad in adinputs:
            # Get the keywords for all sections to be trimmed.
            datasec_kw = ad._keyword_for('data_section')
            darksec_kw = ad._keyword_for('dark_section')

            all_datasecs = ad.data_section()

            for ext, datasec in zip(ad, all_datasecs):
                ext.reset(ext.nddata[:, :, datasec.y1:datasec.y2, datasec.x1:datasec.x2])

                # Update the data section keywords in the header
                sections, new_sections = gt.map_data_sections_to_trimmed_data(datasec)
                xshift = sections[0].x1 - new_sections[0].x1
                yshift = sections[0].y1 - new_sections[0].y1
                ampsec_list = gt.map_array_sections(ext)
                for amp, sec in zip(range(1,len(ampsec_list)+1), ampsec_list):
                    new_x1 = sec.x1 - xshift + 1
                    new_x2 = sec.x2 - xshift
                    new_y1 = sec.y1 - yshift + 1
                    new_y2 = sec.y2 - yshift
                    newDataSecStr = f'[{new_x1}:{new_x2},{new_y1}:{new_y2}]'
                    ext.hdr[f'{datasec_kw}{amp}'] = newDataSecStr

                # Remove the reference pixel and dark section keywords from the headers
                del ext.hdr[f'{darksec_kw}*']
                del ext.hdr[f'REFSCT*']
                del ext.hdr[f'REFSCB*']
                del ext.hdr[f'REFSCS*']

            # Update the filename.
            ad.update_filename(suffix=sfx, strip=True)
        
        return adinputs

    def _modify_filename_components(self, filename, n_int=None, suffix=None):
        """
        Insert a running integration number into the filename, for processing
        integrations that are temporarily split into separate ad instances,
        and/or replace the suffix(es).
        """
        root, fext = os.path.splitext(filename)
        components = root.split("_", maxsplit=1)
        base = components.pop(0)
        suff = ((components[0] if components else '') if suffix is None
                else suffix.lstrip('_'))
        suff = '_' + suff if suff else ''
        return base + ('' if n_int is None else f'-{n_int}') + suff + fext

    def _smoothFFT(self, data, delt, first_deriv=False, second_deriv=False):
        """
        Optimal smoothing algorithm

        Smoothing algorithm to perform optimal filtering of the vertical
        reference pixels to reduce 1/f noise (horizontal stripes), based on the
        Kosarev and Pantos algorithm. This assumes that the data to be 
        filtered/smoothed has been sampled evenly.

        If first_deriv is set, then returns two results.
        If second_deriv is set, then returns three results.

        This function was copied, without modification, from pyNRC. The code
        was accessed between October and December 2020.

        pyNRC:
            https://github.com/JarronL/pynrc
        This code specifically:
            https://github.com/JarronL/pynrc/blob/master/pynrc/reduce/ref_pixels.py

        License
        -------
        MIT License

        Copyright (c) 2018, Jarron Leisenring

        Permission is hereby granted, free of charge, to any person obtaining a copy
        of this software and associated documentation files (the "Software"), to deal
        in the Software without restriction, including without limitation the rights
        to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
        copies of the Software, and to permit persons to whom the Software is
        furnished to do so, subject to the following conditions:

        The above copyright notice and this permission notice shall be included in all
        copies or substantial portions of the Software.

        THIS SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
        IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
        FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
        AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
        LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
        OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
        SOFTWARE.

        Parameters
        ----------
        data : ndarray
            Signal to be filtered.
        delt : float
            Delta time between samples.
        first_deriv : bool
            Return the first derivative.
        second_deriv : bool
            Return the second derivative (along with the first).
        """
        Dat = data.flatten()
        N = Dat.size
        Pi2 = 2 * np.pi
        OMEGA = Pi2 / (N * delt)
        X = np.arange(N) * delt

        ##------------------------------------------------
        ## Center and Baselinefit of the data
        ##------------------------------------------------
        Dat_m = Dat - np.mean(Dat)
        SLOPE = (Dat_m[-1] - Dat_m[0]) / (N-2)
        Dat_b = Dat_m - Dat_m[0] - SLOPE * X / delt

        ##------------------------------------------------
        ## Compute fft- / power- spectrum
        ##------------------------------------------------
        Dat_F = np.fft.rfft(Dat_b)
        Dat_P = np.abs(Dat_F)**2

        ##------------------------------------------------
        ## Noise spectrum from 'half' to 'full'
        ## Mind: half means N/4, full means N/2
        ##------------------------------------------------
        i1 = int((N-1) / 4)
        i2 = int((N-1) / 2) + 1
        Sigma = np.sum(Dat_P[i1:i2])
        Noise = Sigma / ((N-1)/2 - (N-1)/4)

        ##------------------------------------------------
        ## Get Filtercoeff. according to Kosarev/Pantos
        ## Find the J0, start search at i=1 (i=0 is the mean)
        ##------------------------------------------------
        J0 = 2
        for i in np.arange(1, int(N/4)+1):
            sig0, sig1, sig2, sig3 = Dat_P[i:i+4]
            if (sig0<Noise) and ((sig1<Noise) or (sig2<Noise) or (sig3<Noise)):
                J0 = i
                break

        ##------------------------------------------------
        ## Compute straight line extrapolation to log(Dat_P)
        ##------------------------------------------------
        ii = np.arange(1, J0+1)
        logvals = np.log(Dat_P[1:J0+1])
        XY = np.sum(ii * logvals)
        XX = np.sum(ii**2)
        S = np.sum(logvals)
        # Find parameters A1, B1
        XM = (2. + J0) / 2
        YM = S / J0
        A1 = (XY - J0*XM*YM) / (XX - J0*XM*XM)
        B1 = YM - A1 * XM

        # Compute J1, the frequency for which straight
        # line extrpolation drops 20dB below noise
        J1 = int(np.ceil((np.log(0.01*Noise) - B1) / A1))
        if J1<J0:
            J1 = J0+1

        ##------------------------------------------------
        ## Compute the Kosarev-Pantos filter windows
        ## Frequency-ranges: 0 -- J0 | J0+1 -- J1 | J1+1 -- N2
        ##------------------------------------------------
        nvals = int((N-1)/2 + 1)
        LOPT = np.zeros_like(Dat_P)
        LOPT[0:J0+1] = Dat_P[0:J0+1] / (Dat_P[0:J0+1] + Noise)
        i_arr = np.arange(J1-J0) + J0+1
        LOPT[J0+1:J1+1] = np.exp(A1*i_arr+B1) / (np.exp(A1*i_arr+B1) + Noise)

        ##--------------------------------------------------------------------
        ## De-noise the Spectrum with the filter
        ## Calculate the first and second derivative (i.e. multiply by iW)
        ##--------------------------------------------------------------------

        # first loop gives smoothed data
        # second loop produces first derivative
        # third loop produces second derivative
        if second_deriv:
            ndiff = 3
        elif first_deriv:
            ndiff = 2
        else:
            ndiff = 1

        for diff in range(ndiff):
            Fltr_Spectrum = np.zeros_like(Dat_P,dtype=complex)
            # make the filter complex
            i1 = 1; n2 = int((N-1)/2)+1; i2 = i1+n2
            FltrCoef = LOPT[i1:].astype(complex)
            # differentiation in frequency domain
            iW = ((np.arange(n2)+i1)*OMEGA*1j)**diff
            # multiply spectrum with filter coefficient
            Fltr_Spectrum[i1:] = Dat_F[i1:] * FltrCoef * iW

            # Fltr_Spectrum[0] values
            # The derivatives of Fltr_Spectrum[0] are 0
            # Mean if diff = 0
            Fltr_Spectrum[0] = 0 if diff>0 else Dat_F[0]

            # Inverse fourier transform back in time domain
            Dat_T = np.fft.irfft(Fltr_Spectrum)

            # This is the smoothed time series (baseline added)
            if diff==0:
                Smoothed_Data = np.real(Dat_T) + Dat[0] + SLOPE * X / delt
            elif diff==1:
                First_Diff = np.real(Dat_T) + SLOPE / delt
            elif diff==2:
                Secnd_Diff = np.real(Dat_T)

        if second_deriv:
            return Smoothed_Data, First_Diff, Secnd_Diff
        elif first_deriv:
            return Smoothed_Data, First_Diff
        else:
            return Smoothed_Data

    def _split_integrations_to_ad(self, adinputs):
        """
        Split integrations stacked along a third/fourth array axis into
        separate AstroData instances with one axis fewer, so they can be fed to
        core primitives as normal (usually 2D) images. Any inputs that have
        already been reduced to 2D are fed through as single "integrations" in
        the same way, for compatibility.

        If any of the adinputs has a SKYTABLE in the expected format (with
        extra columns for the object & sky integration numbers), it will also
        be split accordingly and propagated to the outputs. Otherwise, only
        tables that are attached to the extensions get propagated.

        Returns
        -------
        tuple : (adoutputs, mapping)
            adoutputs is a list of AstroData instances, derived by splitting
            adinputs into a separate instance per integration (conserving the
            number of extensions, which is normally 1).

            mapping is a dict of { out_filename : (in_ad, n_integ) }
            that maps the filenames of adoutputs to the original adinputs and
            integration numbers, to avoid having to parse these components from
            the modified filenames afterwards.

        """
        adoutputs = []
        mapping = {}

        for ad in adinputs:

            if not ad.filename:
                # Filenames could be made optional, but then we won't be able
                # to generate a useful mapping; do that only if it's needed.
                raise ValueError('AstroData instance has no filename')

            try:
                sky_table = ad.SKYTABLE
            except AttributeError:
                sky_table = None

            n_ints = max(ext.shape[0] if len(ext.shape) > 2 else 1
                         for ext in ad)
            split_ads = n_ints * [None]

            for n_ext, ext in enumerate(ad):

                dt = ad.ut_datetime()  # start of first integration
                tdelt = timedelta(seconds=ext.hdr.get('INTTOTT') or 0.)
                ndd = ext.nddata if len(ext.shape) > 2 else [ext.nddata]

                for n_int, plane in enumerate(ndd):
                    if n_ext == 0:
                        new_ad = split_ads[n_int] = astrodata.create(ad.phu)
                        new_ad.filename = self._modify_filename_components(
                            ad.filename, n_int=n_int
                        )
                        # Could also update MJD-OBS & EXPSTART here, but that
                        # would require AstroPy instead of plain datetime and
                        # currently isn't necessary:
                        new_ad.phu['DATE-OBS'] = dt.date().isoformat()
                        new_ad.phu['TIME-OBS'] = dt.time().isoformat()
                        ## This probably isn't needed, since input HISTORY gets
                        ## updated at the end by the decorator irrespective.
                        # for attr in ad.tables:
                        #     if attr != 'SKYTABLE':
                        #         setattr(new_ad, attr, getattr(ad, attr))
                        if sky_table and 'OBJINT' in sky_table.columns:
                            rows = sky_table[sky_table['OBJINT'] == n_int]
                            new_ad.SKYTABLE = Table(
                                names=('SKYNAME',),
                                dtype=(str,),
                                data = [[
                                    self._modify_filename_components(
                                        row['SKYNAME'], n_int=row['SKYINT']
                                    ) for row in rows
                                ]]
                            )
                        adoutputs.append(new_ad)
                        mapping[new_ad.filename] = (ad, n_int)
                    split_ads[n_int].append(plane, header=ext.hdr)
                    dt += tdelt

        return adoutputs, mapping


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
