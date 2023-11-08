#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_spect_nearIR.py
# ------------------------------------------------------------------------------

import astrodata
from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt

from geminidr.core import NearIR
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_nearIR

from recipe_system.utils.decorators import parameter_override

import numpy as np
import numpy.linalg as la
import sys
from copy import deepcopy
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioNearIR(Scorpio, NearIR):
    """
    This class contains primitives that applies to all Scorpio near-IR data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'NIR'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioNearIR, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR)

    def calculateSignalByRegression(self, adinputs=None, **params):
        """
        Iteratively fit a slope, intercept, and cosmic rays to a ramp.

        This function fits a ramp, possibly with discontinuities (cosmic-ray
        hits), to a 3-D data "cube" with shape (number of groups, number of
        pixels high, number of pixels wide).  The fit will be done multiple
        times, with the previous fit being used to assign weights (via the
        covariance matrix) for the current fit.  The iterations stop either
        when the maximum number of iterations has been reached or when the
        maximum difference between the previous fit and the current fit is
        below a cutoff.  This function calls compute_slope and evaluate_fit.
        
        compute_slope creates arrays for the slope, intercept, and cosmic-ray
        amplitudes (the arrays that will be returned by determine_slope).  Then
        it loops over the number of cosmic rays, from 0 to max_num_cr
        inclusive.  Within this loop, compute_slope copies to temporary arrays
        the ramp data for all the pixels that have the current number of cosmic
        ray hits, calls gls_fit to compute the fit, then copies the results
        of the fit (slope, etc.) to the output arrays for just those pixels.
        The input to gls_fit is ramp data for a subset of pixels (nz in number)
        that all have the same number of cosmic-ray hits.  gls_fit solves
        matrix equations (one for each of the nz pixels) of the form:
            
            y = x * p
        
        where y is a column vector containing the observed data values in
        electrons for each group (the length of y is ngroups, the number of
        groups); x is a matrix with ngroups rows and 2 + num_cr columns, where
        num_cr is the number of cosmic rays being included in the fit; and p
        is the solution, a column vector containing the intercept, slope, and
        the amplitude of each of the num_cr cosmic rays.  The first column of
        x is all ones, for fitting to the intercept.  The second column of x
        is the time (seconds) at the beginning of each group.  The remaining
        num_cr columns (if num_cr > 0) are Heaviside functions, 0 for the
        first rows and 1 for all rows at and following the group containing a
        cosmic-ray hit (each row corresponds to a group).  There will be one
        such column for each cosmic ray, so that the cosmic rays will be fit
        independently of each other.  Whether a cosmic ray hit the detector
        during a particular group was determined by a previous step, and the
        affected groups are flagged in the group data quality array.  In order
        to account for the variance of each observed value and the covariance
        between them (since they're measurements along a ramp), the solution
        is computed in this form (the @ sign represents matrix multiplication):
            
            (xT @ C^-1 @ x)^-1 @ [xT @ C^-1 @ y]
        
        where C is the ngroups by ngroups covariance matrix, ^-1 means matrix
        inverse, and xT is the transpose of x.
        
        Summary of the notation:
        
        data_sect is 3-D, (ngroups, ny, nx); this is the ramp of science data.
        cr_flagged is 3-D, (ngroups, ny, nx); 1 indicates a cosmic ray, e.g.:
            cr_flagged = np.where(np.bitwise_and(gdq_sect, jump_flag), 1, 0)
        cr_flagged_2d is 2-D, (ngroups, nz); this gives the location within
            the ramp of each cosmic ray, for the subset of pixels (nz of them)
            that have a total of num_cr cosmic ray hits at each pixel.  This
            is passed to gls_fit(), which fits a slope to the ramp.
        
        ramp_data has shape (ngroups, nz); this will be a ramp with a 1-D
        array of pixels copied out of data_sect.  The pixels will be those
        that have a particular number of cosmic-ray hits, somewhere within
        the ramp.
        
        Sum cr_flagged over groups to get an (ny, nx) image of the number of
        cosmic rays (i.e. accumulated over the ramp) in each pixel.
        sum_flagged = cr_flagged.sum(axis=0, ...)
        sum_flagged is used to extract the nz pixels from (ny, nx) that have a
        specified number of cosmic ray hits, e.g.:
            for num_cr in range(max_num_cr + 1):
                ncr_mask = (sum_flagged == num_cr)
                nz = ncr_mask.sum(dtype=np.int32)
                for k in range(ngroups):
                    ramp_data[k] = data_sect[k][ncr_mask]
                    cr_flagged_2d[k] = cr_flagged[k][ncr_mask]
        
        gls_fit is called for the subset of pixels (nz of them) that have
        num_cr cosmic ray hits within the ramp, the same number for every
        pixel.

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]
        
        slope_diff_cutoff = 1.e-5

        # Lower and upper limits of the number of iterations that will be done
        # by determineSlope.
        min_iter = 1
        max_iter = 1

        adoutputs = []
        for ad in adinputs:
            new_ad = astrodata.create(ad.phu)

            for e, ext in enumerate(ad):
                ints, groups, rows, cols = ext.data.shape
                new_ext_data = np.empty((ints, rows, cols), dtype=ext.data.dtype)
                new_ext_mask = np.empty((ints, rows, cols), dtype=ext.mask.dtype)
                new_ext_variance = np.empty((ints, rows, cols), dtype=ext.variance.dtype)
                new_ad.append(new_ext_data, header=ext.hdr)

                for intn in range(len(ext.data)):
                    # Get data from this extension.
                    input_var_sect_2 = ext.variance[intn] ** 2
                    gdq_sect = ext.mask[intn]
                    readnoise_sect = gt.array_from_descriptor_value(ext, 'read_noise')[0][0]
                    gain_sect = gt.array_from_descriptor_value(ext, 'gain')[0][0]
                    frame_time = 0.38 #ext.hdr['INTTIME']
                    group_time = frame_time * (ext.hdr['UTRFRAME'] + ext.hdr['UTRSKIP'])
                    nframes = ext.hdr['UTRFRAME']
                    max_num_cr = self._get_max_num_cr(gdq_sect, DQ.cosmic_ray)

                    data = deepcopy(ext.data[intn])

                    # These will be updated in the loop.
                    prev_slope_sect = (data[1, :, :] - data[0, :, :]) / group_time
                    prev_fit = deepcopy(data)

                    iter = 0
                    done = False
                    while not done:
                        (intercept_sect, int_var_sect, slope_sect, slope_var_sect,
                         cr_sect, cr_var_sect) = \
                            self._compute_slope(data, input_var_sect_2, gdq_sect,
                                                readnoise_sect, gain_sect,
                                                prev_fit, prev_slope_sect,
                                                frame_time, group_time, nframes,
                                                max_num_cr, DQ.saturated, DQ.cosmic_ray)
                        iter += 1
                        if iter >= max_iter:
                            done = True
                        else:
                            # If there are pixels with zero or negative variance,
                            # ignore them when taking the difference between the
                            # slopes computed in the current and previous
                            # iterations.
                            slope_diff = np.where(slope_var_sect > 0., prev_slope_sect - slope_sect, 0.)
                            max_slope_diff = np.abs(slope_diff).max()

                            if iter >= min_iter and max_slope_diff < slope_diff_cutoff:
                                done = True

                            current_fit = self._evaluate_fit(intercept_sect, slope_sect, cr_sect,
                                                       frame_time, group_time, gdq_sect, DQ.cosmic_ray)

                            prev_fit = self._positive_fit(current_fit)    # use for next iteration
                            del current_fit

                            prev_slope_sect = slope_sect.copy()

                    # Create a basis for the 2D data quality array
                    temp_dq = np.zeros((gdq_sect.shape[1], gdq_sect.shape[2]), dtype=DQ.datatype)

                    # Replace zero or negative variances with this:
                    LARGE_VARIANCE = 1.e8

                    v_mask = (slope_var_sect <= 0.)
                    if v_mask.any():
                        # Replace negative or zero variances with a large value.
                        slope_var_sect[v_mask] = LARGE_VARIANCE
                        # Also set a flag in the pixel dq array.
                        temp_dq[v_mask] = DQ.bad_pixel
                        del v_mask

                    # If a pixel was flagged (by an earlier step) as saturated in
                    # the first group, flag the pixel as bad.
                    s_mask = (gdq_sect[0] == DQ.saturated)
                    if s_mask.any():
                        temp_dq[s_mask] = DQ.bad_pixel

                    # Compress the DQ array 3D -> 2D
                    pixel_dq = self._dq_compress_sect(gdq_sect, temp_dq)

                    new_ext_data[intn] = slope_sect
                    new_ext_mask[intn] = pixel_dq
                    new_ext_variance[intn] = np.sqrt(slope_var_sect)

                new_ad[e].reset(new_ext_data, mask=new_ext_mask, variance=new_ext_variance)

            # Update the filename.
            ad.update_filename(suffix=sfx, strip=True)
            new_ad.filename = ad.filename
            adoutputs.append(new_ad)

        return adoutputs

    def flagCosmicRaysFromNDRs(self, adinputs=None, **params):
        """
        Two-Point difference method for finding outliers in a 3-D data array.

        The scheme used in this variation of the method uses numpy array methods
        to compute first differences and find the max outlier in each pixel
        while still working in the full 3-D data array. This makes detection of
        the first outlier very fast. We then iterate pixel-by-pixel over only
        those pixels that are already known to contain an outlier, to look for
        any additional outliers and set the appropriate DQ mask for all outliers
        in the pixel.

        This function was copied and modified from stcal. The code was accessed
        between June and July 2021.

        stcal:
            https://github.com/spacetelescope/stcal

        The two-point difference method used in this function is based on the
        method by Anderson & Gordon, 2011: 
            https://iopscience.iop.org/article/10.1086/662593

        TODO:
        - add rej_threshold input
    
        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]

        # TODO - add input to the primitive to set this value
        # This value yields ~0.04% pixels whose largest non-saturated ratio is above the threshold, using test data from PhoSim
        rej_threshold = 50

        for ad in adinputs:
            for e, ext in enumerate(ad):
                for intn in range(len(ext.data)):
                    # Make a copy of the data so we don't make changes to it
                    data = deepcopy(ext.data[intn])

                    # Get the number of frames averaged per group
                    nframes = ext.hdr["UTRFRAME"]

                    # Get data characteristics
                    ngroups, nrows, ncols = data.shape
                    ndiffs = ngroups - 1

                    # Get the read noise array and square it, for use later
                    read_noise = gt.array_from_descriptor_value(ext, "read_noise")[0][0]
                    read_noise_2 = read_noise ** 2

                    # Set saturated values and pixels flagged as bad pixels in the
                    # input data array to NaN, so they don't get used in any of the
                    # subsequent calculations.
                    data[np.where(np.bitwise_and(ext.mask[intn], DQ.saturated))] = np.nan
                    data[np.where(np.bitwise_and(ext.mask[intn], DQ.bad_pixel))] = np.nan

                    # Compute first differences of adjacent groups up the ramp.
                    # Note: roll the ngroups axis of data array to the end, to make
                    # memory access to the values for a given pixel faster.
                    # New form of the array has dimensions [nrows, ncols, ngroups].
                    first_diffs = np.diff(np.rollaxis(data, axis=0, start=3), axis=2)
                    positive_first_diffs = np.abs(first_diffs)

                    # sat_groups is a 3D array that is true when the group is
                    # saturated.
                    sat_groups = np.isnan(positive_first_diffs)

                    # number_sat_groups is a 2D array with the count of saturated
                    # groups for each pixel.
                    number_sat_groups = sat_groups.sum(axis=2)

                    # Make all the first diffs for saturated groups be equal to
                    # 100,000 to put them above the good values in the sorted index.
                    first_diffs[np.isnan(first_diffs)] = 100000.

                    # Here we sort the 3D array along the last axis, which is the
                    # group axis. np.argsort returns a 3D array with the last axis
                    # containing the indices that would yield the groups/diffs in
                    # order.
                    sort_index = np.argsort(positive_first_diffs)

                    # median_diffs is a 2D array with the clipped median of each pixel
                    median_diffs = self._get_clipped_median(ndiffs, number_sat_groups, first_diffs, sort_index)

                    # Compute uncertainties as the quadrature sum of the poisson
                    # noise in the first difference signal and read noise. Because
                    # the first differences can be biased by CRs/jumps, we use the
                    # median signal for computing the poisson noise. Here we lower
                    # the read noise by the square root of number of frames in
                    # the group. Sigma is a 3D array.
                    poisson_noise = np.sqrt(np.abs(median_diffs))
                    poisson_noise_2 = poisson_noise ** 2
                    sigma = np.sqrt(poisson_noise_2 + read_noise_2 / nframes)

                    # Reset sigma to exclude pixels with both readnoise and signal=0
                    sigma_0_pixels = np.where(sigma == 0.)
                    if len(sigma_0_pixels[0] > 0):
                        #log.debug(f'Twopt found {len(sigma_0_pixels[0])} pixels with sigma=0')
                        #log.debug('which will be reset so that no jump will be detected.')
                        huge_num = np.finfo(np.float32).max
                        sigma[sigma_0_pixels] = huge_num

                    # Compute distance of each sample from the median in units of
                    # sigma; note that the use of "abs" means we'll detect positive
                    # and negative outliers. ratio is a 3D array with the units of
                    # sigma deviation of the difference from the median.
                    ratio = np.abs(first_diffs - median_diffs[:, :, np.newaxis]) / sigma[:, :, np.newaxis]

                    # Get the group index for each pixel of the largest
                    # non-saturated group, assuming the indices are sorted. 2 is
                    # subtracted from ngroups because we are using differences
                    # and there is one less difference than the number of groups.
                    # This is a 2D array.
                    max_value_index = ngroups - 2 - number_sat_groups

                    # Extract from the sorted group indices the index of the largest
                    # non-saturated group.
                    row, col = np.where(number_sat_groups >= 0)
                    max_index1d = sort_index[row, col, max_value_index[row, col]]
                    max_index1 = np.reshape(max_index1d, (nrows, ncols)) # reshape to a 2-D array

                    # Get the row and column indices of pixels whose largest
                    # non-saturated ratio is above the threshold.
                    r, c = np.indices(max_index1.shape)
                    row1, col1 = np.where(ratio[r, c, max_index1] > rej_threshold)
                    #log.info(f'From highest outlier Two-point found {len(row1)} pixels with at least one CR')

                    # Loop over all pixels that we found the first CR in
                    number_pixels_with_cr = len(row1)
                    for j in range(number_pixels_with_cr):
                        # Extract the first diffs for this pixel with at least one
                        # CR, yielding a 1-D array.
                        pixel_masked_diffs = first_diffs[row1[j], col1[j]]

                        # Get the scalar readnoise^2 and number of saturated groups
                        # for this pixel
                        pixel_rn2 = read_noise_2[row1[j], col1[j]]
                        pixel_sat_groups = number_sat_groups[row1[j], col1[j]]

                        # Create a CR mask and set 1st CR to be found
                        # cr_mask=0 designates a CR
                        pixel_cr_mask = np.ones(pixel_masked_diffs.shape, dtype=bool)
                        number_CRs_found = 1
                        pixel_sorted_index = sort_index[row1[j], col1[j], :]
                        pixel_cr_mask[pixel_sorted_index[ndiffs - pixel_sat_groups - 1]] = 0 # setting largest diff to be a CR
                        new_CR_found = True

                        # Loop and see if there is more than one CR, setting the
                        # mask as you go
                        while new_CR_found and ((ndiffs - number_CRs_found - pixel_sat_groups) > 1):
                            new_CR_found = False
                            # For this pixel get a new median difference excluding
                            # the number of CRs found and the number of saturated
                            # groups
                            pixel_med_diff = self._get_clipped_median(ndiffs, number_CRs_found + pixel_sat_groups,
                                                                      pixel_masked_diffs, pixel_sorted_index)

                            # Recalculate the noise and ratio for this pixel now
                            # that we have rejected a CR
                            pixel_poisson_noise = np.sqrt(np.abs(pixel_med_diff))
                            pixel_poisson_noise_2 = pixel_poisson_noise ** 2
                            pixel_sigma = np.sqrt(pixel_poisson_noise_2 + pixel_rn2 / nframes)
                            pixel_ratio = np.abs(pixel_masked_diffs - pixel_med_diff) / pixel_sigma

                            # Check if largest remaining difference is above threshold
                            if pixel_ratio[pixel_sorted_index[ndiffs - number_CRs_found - pixel_sat_groups - 1]] > rej_threshold:
                                new_CR_found = True
                                pixel_cr_mask[pixel_sorted_index[ndiffs - number_CRs_found - pixel_sat_groups - 1]] = 0
                                number_CRs_found += 1

                        # Found all CRs for this pixel. Set CR flags in input DQ array for this pixel
                        ext.mask[intn, 1:, row1[j], col1[j]] = \
                            np.bitwise_or(ext.mask[intn, 1:, row1[j], col1[j]], DQ.cosmic_ray * np.invert(pixel_cr_mask))

            # Update the filename.
            ad.update_filename(suffix=sfx, strip=True)

        return adinputs

    def referencePixelsCorrect(self, adinputs=None, **params):
        adinputs = self.subtractReferencePixels(adinputs, 
                    **self._inherit_params(params, "subtractReferencePixels"))
        adinputs = self.trimReferencePixels(adinputs, suffix=params["suffix"])
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

    def _compute_slope(self, data_sect, input_var_sect, gdq_sect,
                       readnoise_sect, gain_sect, prev_fit, prev_slope_sect,
                       frame_time, group_time, nframes, max_num_cr,
                       saturated_flag, jump_flag):
        """
        Set up the call to fit a slope to ramp data.

        This loops over the number of cosmic rays (jumps). That is, all the
        ramps with no cosmic rays are processed first, then all the ramps with
        one cosmic ray, then with two, etc.

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        data_sect : 3-D ndarray; shape (ngroups, ny, nx)
            The ramp data for one extension of an ad object. This may be a
            subarray in detector coordinates, but covering all groups.
        
        input_var_sect: 3-D ndarray; shape (ngroups, ny, nx)
            The square of the input variance array, matching data_sect.
        
        gdq_sect : 3-D ndarray; shape (ngrops, ny, nx)
            The group data quality array. This may be a subarray, matching
            data_sect.
        
        readnoise_sect : 2-D ndarray; shape (ny, nx)
            The read noise in electrons at each detector pixel (i.e. not a
            ramp). This may be a subarray, similar to data_sect.
        
        gain_sect : 2-D ndarray or None; shape (ny, nx)
            The gain in electrons per DN at each detector pixel (i.e. not a
            ramp). This may be a subarray, matching readnoise_sect. If gain_sect
            is None, a value of 1 will be assumed.
        
        prev_fit : 3-D ndarray; shape (ngroups, ny, nx)
            The previous fit (intercept, slope, cosmic-ray amplitudes)
            evaluated for each pixel in the subarray.  data_sect itself may be
            used for the first iteration.
        
        prev_slope_sect : 2-D ndarray; shape (ny, nx)
            An estimate (e.g. from a previous iteration) of the slope at each
            pixel, in electrons per second.  This may be a subarray, similar to
            data_sect.
        
        frame_time : float
            The time to read one frame, in seconds.
        
        group_time : float
            Time increment between groups, in seconds.
        
        nframes : int
            Number of frames that were averaged together to make a group.
            This value does not include the number (if any) of skipped frames.
        
        max_num_cr : non-negative int
            The maximum number of cosmic rays that should be handled.
        
        saturated_flag : int
            DQ_definitions['saturated']
        
        jump_flag : int
            DQ_definitions['cosmic_ray']

        Returns
        -------
        tuple : (intercept_sect, int_var_sect, slope_sect, slope_var_sect,
                 cr_sect, cr_var_sect)
            intercept_sect is a 2-D ndarray (ny, nx), the intercept of the ramp
            at each pixel of data_sect.

            int_var_sect is a 2-D ndarray (ny, nx), the variance of the
            intercept at each pixel of data_sect.

            slope_sect is a 2-D ndarray (ny, nx), the ramp slope at each pixel
            of data_sect.

            slope_var_sect is a 2-D ndarray (ny, nx), the variance of the slope
            at each pixel of data_sect.

            cr_sect is a 3-D ndarray (ny, nx, cr_dimen), the amplitude of each
            cosmic ray at each pixel of data_sect. cr_dimen is max_num_cr or 1,
            whichever is larger.

            cr_var_sect is a 3-D ndarray (ny, nx, cr_dimen), the variance of
            each cosmic ray amplitude.
        """
        cr_flagged = np.empty(data_sect.shape, dtype=np.uint8)
        cr_flagged[:] = np.where(np.bitwise_and(gdq_sect, jump_flag), 1, 0)

        # If a pixel is flagged as a jump in the first group, we can't fit to
        # the ramp, because a matrix that we need to invert would be singular.
        # If there's only one group, we can't fit a ramp to it anyway, so
        # at this point we wouldn't need to be concerned about a jump. If
        # there is more than one group, just ignore any jump in the first group.
        if data_sect.shape[0] > 1:
            cr_flagged[0, :, :] = 0

        # Sum over groups to get an (ny, nx) image of the number of cosmic
        # rays in each pixel, accumulated over the ramp.
        sum_flagged = cr_flagged.sum(axis=0, dtype=np.int32)

        # If a pixel is flagged as saturated in the first or second group, we
        # don't want to even attempt to fit a slope to the ramp for that pixel.
        # Handle this case by setting the corresponding pixel in sum_flagged to
        # a negative number. The test 'ncr_mask = (sum_flagged == num_cr)'
        # will therefore never match, since num_cr is zero or larger, and the
        # pixel will not be included in any ncr_mask.
        mask1 = (gdq_sect[0, :, :] == saturated_flag)
        sum_flagged[mask1] = -1

        # one_group_mask flags pixels that are not saturated in the first
        # group but are saturated in the second group (if there is a second
        # group). For these pixels, we will assign a value to the slope
        # image by just dividing the value in the first group by group_time.
        if len(gdq_sect) > 1:
            mask2 = (gdq_sect[1, :, :] == saturated_flag)
            sum_flagged[mask2] = -1
            one_group_mask = np.bitwise_and(mask2, np.bitwise_not(mask1))
            del mask2
        else:
            one_group_mask = np.bitwise_not(mask1)
        del mask1

        # Set elements of this array to a huge value if the corresponding pixels
        # are saturated. This is not a flag, it's a value to be added to the
        # diagonal of the covariance matrix. 1.e20
        saturated = np.empty(data_sect.shape, dtype=np.float64)
        saturated[:] = np.where(np.bitwise_and(gdq_sect, saturated_flag), 1.e20, 0)

        # Create arrays to be populated and then returned.
        shape = data_sect.shape
        dtype = data_sect.dtype
        # Lower limit of one, in case there are no cosmic rays at all.
        cr_dimen = max(1, max_num_cr)
        intercept_sect = np.zeros((shape[1], shape[2]), dtype=dtype)
        slope_sect = np.zeros((shape[1], shape[2]), dtype=dtype)
        cr_sect = np.zeros((shape[1], shape[2], cr_dimen), dtype=dtype)
        int_var_sect = np.zeros((shape[1], shape[2]), dtype=dtype)
        slope_var_sect = np.zeros((shape[1], shape[2]), dtype=dtype)
        cr_var_sect = np.zeros((shape[1], shape[2], cr_dimen), dtype=dtype)

        # This takes care of the case that there's only one group, as well as
        # pixels that are saturated in the second but not the first group.
        if one_group_mask.any():
            slope_sect[one_group_mask] = data_sect[0, one_group_mask] / group_time
        del one_group_mask

        # Fit slopes for all pixels that have no cosmic ray hits anywhere in
        # the ramp, then fit slopes with one CR hit, then with two, etc.
        for num_cr in range(max_num_cr + 1):
            ngroups = len(data_sect)
            ncr_mask = (sum_flagged == num_cr)

            # Number of detector pixels flagged with num_cr CRs within the ramp.
            nz = ncr_mask.sum(dtype=np.int32)
            if nz <= 0:
                continue

            # ramp_data will be a ramp with a 1-D array of pixels copied out
            # of data_sect.
            ramp_data = np.empty((ngroups, nz), dtype=data_sect.dtype)
            input_var_data = np.empty((ngroups, nz), dtype=data_sect.dtype)
            prev_fit_data = np.empty((ngroups, nz), dtype=prev_fit.dtype)
            prev_slope_data = np.empty(nz, dtype=prev_slope_sect.dtype)
            readnoise = np.empty(nz, dtype=readnoise_sect.dtype)

            prev_slope_data[:] = prev_slope_sect[ncr_mask]
            readnoise[:] = readnoise_sect[ncr_mask]
            
            if gain_sect is None:
                gain = None
            else:
                gain = np.empty(nz, dtype=gain_sect.dtype)
                gain[:] = gain_sect[ncr_mask]
            
            cr_flagged_2d = np.empty((ngroups, nz), dtype=cr_flagged.dtype)
            saturated_data = np.empty((ngroups, nz), dtype=prev_fit.dtype)
            
            for k in range(ngroups):
                ramp_data[k] = data_sect[k][ncr_mask]
                input_var_data[k] = input_var_sect[k][ncr_mask]
                prev_fit_data[k] = prev_fit[k][ncr_mask]
                cr_flagged_2d[k] = cr_flagged[k][ncr_mask]
                # This is for clobbering saturated pixels.
                saturated_data[k] = saturated[k][ncr_mask]

            (result, variances) = \
                self._gls_fit(ramp_data, prev_fit_data, prev_slope_data,
                              readnoise, gain, frame_time, group_time,
                              nframes, num_cr, cr_flagged_2d, saturated_data)

            # Copy the intercept, slope, and cosmic-ray amplitudes and their
            # variances to the arrays to be returned.
            # ncr_mask is a mask array that is True for each pixel that has the
            # current number (num_cr) of cosmic rays. Thus, the output arrays
            # are being populated here in sets, a different set of pixels
            # with each iteration of this loop.
            intercept_sect[ncr_mask] = result[:, 0].copy()
            int_var_sect[ncr_mask] = variances[:, 0].copy()
            slope_sect[ncr_mask] = result[:, 1].copy()
            slope_var_sect[ncr_mask] = variances[:, 1].copy()
            
            # In this loop, i is just an index. cr_sect is populated for
            # number of cosmic rays = 1 to num_cr, inclusive.
            for i in range(num_cr):
                cr_sect[ncr_mask, i] = result[:, 2 + i].copy()
                cr_var_sect[ncr_mask, i] = variances[:, 2 + i].copy()

        return (intercept_sect, int_var_sect, slope_sect, slope_var_sect, cr_sect, cr_var_sect)

    def _dq_compress_sect(self, gdq_sect, pixel_dq):
        """
        Compress a 3D data quality cube into a 2D data quality array.

        Parameters
        ----------
        gdq_sect : 3-D ndarray
            The 3D data quality cube for a data section.

        Returns
        -------
        pixel_dq : 2-D ndarray
            The 2D data quality plane for a cube compressed to a 2D array.
        """
        """
        loc_ramp = np.bitwise_and(gdq_sect, DQ.bad_pixel)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.bad_pixel)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.non_linear)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.non_linear)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.saturated)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.saturated)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.cosmic_ray)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.cosmic_ray)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.no_data)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.no_data)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.overlap)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.overlap)

        loc_ramp = np.bitwise_and(gdq_sect, DQ.unilluminated)
        loc_image = np.where(loc_ramp.sum(axis=0) > 0)
        pixel_dq[loc_image] = np.bitwise_or(pixel_dq[loc_image], DQ.unilluminated)
        """

        pixel_dq = np.bitwise_or(np.bitwise_or.reduce(gdq_sect, axis=0), pixel_dq)

        return pixel_dq

    def _evaluate_fit(self, intercept_sect, slope_sect, cr_sect,
                      frame_time, group_time, gdq_sect, jump_flag):
        """
        Evaluate the fit (intercept, slope, cosmic-ray amplitudes).

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        intercept_sect : 2-D ndarray
            The intercept of the ramp at each pixel of data_sect (one of the
            arguments to determine_slope).

        slope_sect : 2-D ndarray
            The ramp slope at each pixel of data_sect.

        cr_sect : 3-D ndarray
            The amplitude of each cosmic ray at each pixel of data_sect.

        frame_time : float
            The time to read one frame, in seconds.

        group_time : float
            Time increment between groups, in seconds.

        gdq_sect : 3-D ndarray; indices:  group, y, x
            The group data quality array.  This may be a subarray, matching
            data_sect.

        jump_flag : int
            The data quality flag indicating a cosmic-ray hit.

        Returns
        -------
        fit_model : 3-D ndarray, shape (ngroups, ny, nx)
            This is the same sahpe as data_sect, and if the fit is good,
            fit_model and data_sect should not differ by much.
        """

        shape_3d = gdq_sect.shape    # the ramp, (ngroups, ny, nx)
        ngroups = gdq_sect.shape[0]

        # This array is also created in the function _compute_slope.
        cr_flagged = np.empty(shape_3d, dtype=np.uint8)
        cr_flagged[:] = np.where(np.bitwise_and(gdq_sect, jump_flag), 1, 0)

        sum_flagged = cr_flagged.sum(axis=0, dtype=np.int32)

        # local_max_num_cr is local to this function. It may be smaller than
        # the max_num_cr that's computed in determine_slope, and it can even be
        # zero.
        local_max_num_cr = sum_flagged.max()
        del sum_flagged

        # The independent variable, in seconds at each image pixel.
        ind_var = np.zeros(shape_3d, dtype=np.float64)
        M = round(group_time / frame_time)
        iv = np.arange(ngroups, dtype=np.float64) * group_time + frame_time * (M + 1.) / 2.
        iv = iv.reshape((ngroups, 1, 1))
        ind_var += iv

        # No cosmic rays yet; these will be accounted for below.
        # ind_var has a different shape (ngroups, ny, nx) from slope_sect and
        # intercept_sect, but their last dimensions are the same.
        fit_model = ind_var * slope_sect + intercept_sect

        # heaviside and cr_flagged have shape (ngroups, ny, nx).
        heaviside = np.zeros(shape_3d, dtype=np.float64)
        cr_cumsum = cr_flagged.cumsum(axis=0, dtype=np.int16)

        # Add an offset for each cosmic ray.
        for n in range(local_max_num_cr):
            heaviside[:] = np.where(cr_cumsum > n, 1., 0.)
            fit_model += (heaviside * cr_sect[:, :, n])

        return fit_model
    
    def _get_clipped_median(self, num_differences, diffs_to_ignore, differences, sorted_index):
        """
        This routine will return the clipped median for the input array or
        pixel. It will ignore the input number of largest differences. At a
        minimum this is at least one plus the number of saturated values, to
        avoid the median being biased by a cosmic ray. As cosmic rays are found,
        the diff_to_ignore will increase.

        This function was copied, without modification from stcal. The code was
        accessed between June and July 2021.

        stcal:
            https://github.com/spacetelescope/stcal

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.
        """

        # Ignore largest value and number of CRs found when finding new median
        # Check to see if this is a 2D array or 1D
        if sorted_index.ndim > 1:
            # Get the index of the median value always excluding the highest
            # value. In addition, decrease the index by 1 for every two
            # diffs_to_ignore, these will be saturated values in this case.
            row, col = np.indices(diffs_to_ignore.shape)
            pixel_med_index = sorted_index[row, col, (num_differences - (diffs_to_ignore[row, col] + 1)) // 2]
            pixel_med_diff = differences[row, col, pixel_med_index]

            # For pixels with an even number of differences the median is the
            # mean of the two central values. So we need to get the value of the
            # other central difference one lower in the sorted index than the
            # one found above.
            even_group_rows, even_group_cols = np.where((num_differences - diffs_to_ignore - 1) % 2 == 0)
            pixel_med_index2 = np.zeros_like(pixel_med_index)
            pixel_med_index2[even_group_rows, even_group_cols] = \
                sorted_index[even_group_rows, even_group_cols, 
                             (num_differences - (diffs_to_ignore[even_group_rows, even_group_cols] + 3)) // 2 ]

            # Average together the two central values
            pixel_med_diff[even_group_rows, even_group_cols] = (
                pixel_med_diff[even_group_rows, even_group_cols] +
                differences[even_group_rows, even_group_cols, pixel_med_index2[even_group_rows, even_group_cols]]) / 2.0

        # The 1-D array case is a lot simpler
        else:
            pixel_med_index = sorted_index[int(((num_differences - 1 - diffs_to_ignore) / 2))]
            pixel_med_diff = differences[pixel_med_index]
            if (num_differences - diffs_to_ignore - 1) % 2 == 0:    # even number of differences
                pixel_med_index2 = sorted_index[int((num_differences - 1 - diffs_to_ignore) / 2) - 1]
                pixel_med_diff = (pixel_med_diff + differences[pixel_med_index2]) / 2.0

        return pixel_med_diff

    def _get_max_num_cr(self, gdq_cube, jump_flag):
        """
        Find the maximum number of cosmic-ray hits in any one pixel.

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        gdq_cube : ndarray
            The group data quality array, 3-D flag
        jump_flag : int
            The data quality flag indicating a cosmic-ray hit.

        Returns
        -------
        max_num_cr : int
            The maximum number of cosmic-ray hits for any pixel.
        """
        cr_flagged = np.empty(gdq_cube.shape, dtype=np.uint8)
        cr_flagged[:] = np.where(np.bitwise_and(gdq_cube, jump_flag), 1, 0)
        max_num_cr = cr_flagged.sum(axis=0, dtype=np.int32).max()

        return max_num_cr

    def _gls_fit(self, ramp_data, prev_fit_data, prev_slope_data,
                 readnoise, gain, frame_time, group_time, nframes,
                 num_cr, cr_flagged_2d, saturated_data):
        """
        Generalized least squares linear fit.

        It is assumed that every input pixel has num_cr cosmic-ray hits
        somewhere within the ramp. This function should be called separately
        for different values of num_cr.

        Notes
        -----
        Currently the noise model is assumed to be a combination of read and
        photon noise alone. Same technique could be used with more complex
        noise models, but then the ramp covariance matrix should be input.

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        ramp_data : 2-D ndarray; indices: group, pixel number
            The ramp data for one of the integrations in an exposure. This may
            be a subset in detector coordiantes, but covering all groups. The
            shape is (ngroups, nz), where ngroups is the length of the ramp,
            and nz is the number of pixels in the current subset.

        prev_fit_data : 2-D ndarray, shape (ngroups, nz)
            The fit to ramp_data, based on applying the values of intercept,
            slope, and cosmic-ray amplitudes that were determined in a previous
            call to gls_fit. This array is only used for setting up the
            covariance matrix.

        prev_slope_data : 1-D ndarray, length nz
            An estimate (ee.g. from a previous iteration) of the slope at each
            pixel, in electrons per second.

        readnoise : 1-D ndarray, length nz
            The read noise in electrons at each detector pixel.

        gain : 1-D ndarray, shape (nz,)
            The analog-to-digital gain (electrons per ADU) at each detector
            pixel.

        frame_time : float
            The time to read one framem in seconds.

        group_time: float
            Time increment between groups, in seconds.

        nframes : int
            Number of frames that were averaged together to make a group.
            This value does not include the number (if any) of skipped frames.

        num_cr : int
            The number of cosmic rays that will be handled. All pixels in the
            current set (ramp_data) are assumed to have this many cosmic ray
            hits somewhere within the ramp.

        cr_flagged_2d: 2-D ndarray, shape (ngroups, nz)
            The values should be 0 or 1; 1 indicates that a cosmic ray was
            detected (by another step) at that point.

        saturated_data : 2-D ndarray, shape (ngroups, nz)
            Normal values are zero; the value will be a huge number for
            saturated pixels. This will be added to the main diagonal of the
            inverse of the weight matrix to greatly reduce the weight for
            saturated pixels.

        Returns
        -------
        tuple : (result2s, variances)
            result2d is a 2-D ndarray; shape (nz, 2 + num_cr)
            The computed values of intercept, slope, and cosmic-ray amplitudes
            (there will be num_cr cosmic-ray amplitudes) for each ofthe nz
            pixels.

            variances is a 2-D ndarray; shape (nz, 2 + num_cr)
            The variance for the intercept, slope, and for the amplitude of
            each cosmic ray that detected.
        """
        M = float(nframes)
        
        ngroups = ramp_data.shape[0]
        nz = ramp_data.shape[1]
        num_cr = int(num_cr)
        
        # x is an array (length nz) of matrices, each of which is the
        # independent variable of a linear equation.  Each such matrix
        # has ngroups rows and 2 + num_cr columns.  The first column is set
        # to 1, for finding the intercept.  The second column is the time at
        # each group, for finding the slope.  The remaining columns (if any),
        # are 0 for all rows prior to a certain point, then 1 for all
        # subsequent rows (i.e. the Heaviside function).  The transition from
        # 0 to 1 is the location of a cosmic ray hit; the first 1 in a column
        # corresponds to the value in cr_flagged_2d being 1.
        x = np.zeros((nz, ngroups, 2 + num_cr), dtype=np.float64)
        x[:, :, 0] = 1.
        x[:, :, 1] = np.arange(ngroups, dtype=np.float64) * group_time + \
            frame_time * (M + 1.) / 2.
        
        if num_cr > 0:
            sum_crs = cr_flagged_2d.cumsum(axis=0)
            for k in range(ngroups):
                s = slice(k, ngroups)
                for n in range(1, num_cr + 1):
                    temp = np.where(np.logical_and(cr_flagged_2d[k] == 1,
                                                   sum_crs[k] == n))
                    if len(temp[0]) > 0:
                        index = (temp[0], s, n + 1)
                        x[index] = 1
            del temp, index
        
        y = np.transpose(ramp_data, (1, 0)).reshape((nz, ngroups, 1))
        
        # ramp_cov is an array of nz matrices, each ngroups x ngroups.
        # each matrix gives the covariance of that pixel's ramp data
        ramp_cov = np.ones((nz, ngroups, ngroups), dtype=np.float64)
        
        # Use the previous fit to the data to populate the covariance matrix,
        # for each of the nz pixels.  prev_fit_data has shape (ngroups, nz),
        # similar to the ramp data, but we want the nz axis to be the first
        # (we're constructing an array of nz matrix equations), so transpose
        # prev_fit_data.
        prev_fit_T = np.transpose(prev_fit_data, (1, 0))
        for k in range(ngroups):
            # Populate the upper right, row by row.
            ramp_cov[:, k, k:ngroups] = prev_fit_T[:, k:k + 1]
            # Populate the lower left, column by column.
            ramp_cov[:, k:ngroups, k] = prev_fit_T[:, k:k + 1]
            # Give saturated pixels a very high high variance (hence a low weight)
            ramp_cov[:, k, k] += saturated_data[k, :]
        del prev_fit_T
        
        # iden is 2-D, but it can broadcast to 4-D.  This is used to add terms to
        # the diagonal of the covariance matrix.
        iden = np.identity(ngroups)
        
        rn3d = readnoise.reshape((nz, 1, 1))
        ramp_cov += (iden * rn3d**2)
        
        # prev_slope_data must be non-negative.
        flags = prev_slope_data < 0.
        prev_slope_data[flags] = 1.
        
        # The resulting fit parameters are
        #  (xT @ ramp_cov^-1 @ x)^-1 @ [xT @ ramp_cov^-1 @ y]
        #  = [y-intercept, slope, cr_amplitude_1, cr_amplitude_2, ...]
        # where @ means matrix multiplication.
        
        # shape of xT is (nz, 2 + num_cr, ngroups)
        xT = np.transpose(x, (0, 2, 1))
        
        # shape of `ramp_invcov` is (nz, ngroups, ngroups)
        iden = iden.reshape((1, ngroups, ngroups))
        ramp_invcov = la.solve(ramp_cov, iden)
        
        del iden
        
        # temp1 = xT @ ramp_invcov
        # np.einsum use is equivalent to matrix multiplication
        # shape of temp1 is (nz, 2 + num_cr, ngroups)
        temp1 = np.einsum('...ij,...jk->...ik', xT, ramp_invcov)
        
        # temp_var = xT @ ramp_invcov @ x
        # shape of temp_var is (nz, 2 + num_cr, 2 + num_cr)
        temp_var = np.einsum('...ij,...jk->...ik', temp1, x)
        
        # `fitparam_cov` is an array of nz covariance matrices.
        # fitparam_cov = (xT @ ramp_invcov @ x)^-1
        # shape of fitparam_covar is (nz, 2 + num_cr, 2 + num_cr)
        I_2 = np.eye(2 + num_cr).reshape((1, 2 + num_cr, 2 + num_cr))
        try:
            # inverse of temp_var
            fitparam_cov = la.solve(temp_var, I_2)
        except la.LinAlgError:
            # find the pixel with the singular matrix
            for z in range(nz):
                try:
                    la.solve(temp_var[z], I_2)
                except la.LinAlgError as msg2:
                    log.warning("singular matrix, z = %d" % z)
                    raise la.LinAlgError(msg2)
        del I_2
        
        # [xT @ ramp_invcov @ y]
        # shape of temp2 is (nz, 2 + num_cr, 1)
        temp2 = np.einsum('...ij,...jk->...ik', temp1, y)
        
        # shape of fitparam is (nz, 2 + num_cr, 1)
        fitparam = np.einsum('...ij,...jk->...ik', fitparam_cov, temp2)
        r_shape = fitparam.shape
        fitparam2d = fitparam.reshape((r_shape[0], r_shape[1]))
        del fitparam
        
        # shape of both result2d and variances is (nz, 2 + num_cr)
        fitparam_uncs = fitparam_cov.diagonal(axis1=1, axis2=2).copy()
        
        return (fitparam2d, fitparam_uncs)

    def _positive_fit(self, current_fit):
        """
        Replace zero and negative values with a positive number.

        Ramp data should be positive, since they are based on counts. The fit to
        a ramp can go negative, however, due e.g. to extrapolation beyond where
        the data are saturated. To avoid negative elements in the covariance
        matrix (populated in part with the fit to the ramp), this function
        replaces zero or negative values in the fit with a positive number.

        License
        -------
        Copyright (C) 2020 Association of Universities for Research in Astronomy (AURA)

        Redistribution and use in source and binary forms, with or without
        modification, are permitted provided that the following conditions are met:

            1. Redistributions of source code must retain the above copyright
              notice, this list of conditions and the following disclaimer.

            2. Redistributions in binary form must reproduce the above
              copyright notice, this list of conditions and the following
              disclaimer in the documentation and/or other materials provided
              with the distribution.

            3. The name of AURA and its representatives may not be used to
              endorse or promote products derived from this software without
              specific prior written permission.

        THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
        WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
        MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
        DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
        INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
        BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
        OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
        ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
        TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
        USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
        DAMAGE.

        Parameters
        ----------
        current_fit : 3-D ndarray, shape (ngroups, ny, nx)
            The fit returned by _evaluate_fit.

        Returns
        -------
        current_fit : 3-D ndarray, shape (ngroups, ny, nx)
            This is the same as the input current_fit, except that zero and
            negative values will have been replaced by a positive value.
        """

        # This is a value to replace zero or negative values in a fit, to make
        # all values of the fit positive and to give low weight where the fit
        # was zero or negative.
        fit_must_be_positive = 1.e10

        return np.where(current_fit <= 0., fit_must_be_positive, current_fit)

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