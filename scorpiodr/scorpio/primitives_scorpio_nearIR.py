#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_spect_nearIR.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt

from geminidr.core import NearIR
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_nearIR

from recipe_system.utils.decorators import parameter_override

import numpy as np
import sys
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

    def detectJumps(self, adclass=None, **params):
        """
        Two-Point difference method for finding outliers in a 3-D data array.

        The scheme used in this variation of the method uses numpy array methods
        to compute first differences and find the max outlier in each pixel
        while still working in the full 3-D data array. This makes detection of
        the first outlier very fast. We then iterate pixel-by-pixel over only
        those pixels that are already known to contain an outlier, to look for
        any additional outliers and set the appropriate DQ mask for all outliers
        in the pixel.

        This function requires at least three good groups in the data ramp.

        This function has been modified from its original location, here:
        https://github.com/spacetelescope/stcal
        Accessed June 2021-

        The two-point difference method used in this function is based on the
        method by Anderson & Gordon, 2011: 
            https://iopscience.iop.org/article/10.1086/662593

        TODO: add AURA lciensing info (https://github.com/spacetelescope/stcal/blob/main/LICENSE)
    
        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        sfx = params["suffix"]

        for ad in adinputs:
            for e, ext in enumerate(ad):
                # Get data characteristics
                ngroups, nrows, ncols = ext.data.shape
                ndiffs = ngroups - 1

                # Get the read noise array and square it, for use later
                read_noise = gt.array_from_descriptor_value(ext, "read_noise")
                read_noise_2 = read_noise ** 2

                # Set saturated values and pixels flagged as bad pixels in the
                # input data array to NaN, so they don't get used in any of the
                # subsequent calculations.
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.saturated))] = np.nan
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.bad_pixel))] = np.nan

                # Compute first differences of adjacent groups up the ramp.
                # Note: roll the ngroups axis of data array to the end, to make
                # memory axxess to the values for a given pixel faster.
                # New form of the array has dimensions [nrows, ncols, ngroups].
                first_diffs = np.diff(np.rollaxis(ext.data, axis=0, start=3), axis=2)
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
                # containing the indices that would yield the groups in order.
                sort_index = np.argsort(positive_first_diffs)

                # median_diffs is a 2D array with the clipped median of each pixel
                median_diffs = self._get_clipped_median_array(ndiffs, number_sat_groups, first_diffs, sort_index)

                # Currently we do not have a keyword for nframes, the number of samples
                # averaged together into each frame. 
                # nframes is the number of samples averaged into a single group.
                #nframes = 

                # Compute uncertainties as the quadrature sum of the poisson
                # noise in the first difference signal and read noise. Because
                # the first differences can be biased by CRs/jumps, we use the
                # median signal for computing the poisson noise. Here we lower
                # the read noise by the square root of number of frames in
                # the group. Sigma is a 3D array.
                sigma = np.sqrt(np.abs(median_diffs) + read_noise_2 / nframes)

                # Reset sigma to exclude pixels with both readnoise and signal=0
                sigma_0_pixels = np.where(sigma == 0.)
                if len(sigma_0_pixels[0] > 0):
                    #log.debug(f'Twopt found {len(sigma_0_pixels[0])} pixels with sigma=0')
                    #log.debug('which will be reset so that no jump will be detected.')
                    huge_num = np.finfo(np.float32).max
                    sigma[sigma_0_pixels] = huge_num

                # Compute distance of each sample from the median in units of
                # sigma; note that the use of "abs" means we'll detect positive
                # and negative outliers.
                # ratio is a 2D array with the units of sigma deviation of the
                # difference from the median.
                ratio = np.abs(first_diffs - median_diffs[:, :, np.newaxis]) / sigma[:, :, np.newaxis]
                ratio3d = np.reshape(ratio, (nrows, ncols, ndiffs))

                # Get the group index for each pixel of the largest
                # non-saturated group, assuming the indices are sorted. 2 is
                # subtracted from ngroups because we are using differences
                # and there is one less difference than the number of groups.
                # This is a 2D array.
                max_value_index = ngroups - 2 - number_sat_groups

        return adinputs






                

    def referencePixelsCorrect(self, adinputs=None, **params):
        adinputs = self.subtractReferencePixels(adinputs, 
                    **self._inherit_params(params, "subtractReferencePixels"))
        adinputs = self.trimReferencePixels(adinputs, suffix=params["suffix"])
        return adinputs

    def subtractReferencePixels(self, adinputs=None, **params):
        """
        Correct a SCORPIO NIR image's noise by using its reference pixels. 

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

        for ad in adinputs:
            # Get the list of section values for all reference pixel sections. 
            refpix_sec = ad.refpix_section()

            for e, ext in enumerate(ad):
                # Convert any saturated or bad pixels to NaNs so they are not
                # used in any calculations.
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.saturated))] == np.nan
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.bad_pixel))] == np.nan

                # Extract the horizontal reference pixel sections.
                rpix_top_sec = refpix_sec['top'][e]
                rpix_bot_sec = refpix_sec['bottom'][e]

                # Get the horizontal reference pixel sections from the data 
                # plane. rpix_top and rpix_bot are lists of 3D arrays 
                # with shape (nframes, nrows, ncols), where each array
                # corresponds to an amplifier.
                rpix_top = [ext.data[:, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_top_sec]
                rpix_bot = [ext.data[:, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_bot_sec]

                # These just make loop indices easier later.
                nframes, nrows, ncols = ext.data.shape
                namps = len(rpix_top)

                # Exctract the vertical reference pixel sections.
                rpix_side_sec = refpix_sec['side'][e]

                # Get the vertical reference pixel data sections from the 
                # data plane. rpix_side is a list of 3D arrays with shape
                # (nframes, nrows, ncols). Note, these arrays do not yet
                # correspond to amplifiers.
                rpix_side = [ext.data[:, sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_side_sec]

                # Restructure the vertical reference pixel sections into 2 
                # sections (left and right). rpix_side_fixed is a list of 
                # 3D arrays with shape (nframes, nrows, ncols). These arrays
                # correspond to amplifiers.
                # This only satisfies full FOV imaging mode.
                rpix_left, rpix_right = [], []
                for sec in range(len(rpix_side)):
                    if sec % 2 == 0:
                        rpix_right.append(rpix_side[sec])
                    else:
                        rpix_left.append(rpix_side[sec])
                rpix_left = np.hstack(rpix_left)
                rpix_right = np.hstack(rpix_right)
                rpix_side_fixed = [rpix_left, rpix_right]

                # Determine the size and shape of one amplifier.
                amp_rows = rpix_side_fixed[0].shape[1] + 8    # +8 for top and bottom refpixels (there will always be 4 on top and 4 on bottom)
                amp_cols = rpix_top_sec[0].x2 - rpix_top_sec[0].x1
                amp_shape = (amp_rows, amp_cols)
                amp_size = amp_rows * amp_cols

                # Create a 2D time array that matches the size and shape of one 
                # amplifier. Multiplly it by 10E-6 (seconds) for pixel readout 
                # rate.
                #amp_shape = (2048, 64) # nrows x ncols
                #amp_size = amp_shape[0] * amp_shape[1]
                time = np.linspace(1, amp_size, amp_size, endpoint=True).reshape((amp_shape[0], -1)) * 10E-6
                time = np.flip(time, 0)

                # Get the list of data section values. We need the indices of 
                # each amplifier in the data section rather than one whole data 
                # section.
                datasec = gt.map_array_sections(ext)

                # Determine the SUPERAVERAGE of the reference pixels. The 
                # superaverage is the average from the top and bottom 
                # reference pixel sections, across all amplifiers, and across 
                # all frames.
                superaverage = np.mean([rpix_top, rpix_bot])

                # Loop over the frames and amplifiers and calculate the average 
                # from the top and bottom reference pixels. This is the 
                # amplifier average. Subtract the superaverage from the 
                # amplifier average to compute the delta. Subtract the delta 
                # from all pixels in this amplifier in this frame.
                for frm in range(nframes):
                    amp_deltas = []
                    for amp in range(namps):
                        ampavg = np.mean([rpix_top[amp][frm], rpix_bot[amp][frm]])
                        delta = ampavg - superaverage
                        amp_deltas.append(delta)

                        rpix_top[amp][frm] -= delta
                        rpix_bot[amp][frm] -= delta

                    # Subtract the delta from the side reference pixel sections.
                    rpix_side_fixed[0][frm] -= amp_deltas[0]     # first amplifier
                    rpix_side_fixed[-1][frm] -= amp_deltas[-1]   # last amplifier

                    # Subtract the delta from the data pixel sections.
                    for amp in range(len(datasec)):
                        ext.data[frm, datasec[amp].y1:datasec[amp].y2, datasec[amp].x1:datasec[amp].x2] -= amp_deltas[amp+1]

                # From the time array, which is the size of the original
                # amplifier, slice out the corresponding horizontal and vertical
                # reference pixels.
                time_top = time[rpix_top_sec[0].y1:rpix_top_sec[0].y2, 0:amp_cols]
                time_bot = time[rpix_bot_sec[0].y1:rpix_bot_sec[0].y2, 0:amp_cols]
                time_side = time[rpix_side_sec[2].y1:-rpix_side_sec[2].y1, rpix_side_sec[0].x1:rpix_side_sec[0].x2]

                # Create lists for the arrays of horizontal and vertical ramp
                # averages.
                ramp_avg_top = []
                ramp_avg_bot = []
                ramp_avg_side = []

                # Loop over the horizontal reference pixels to calculate their
                # ramp averages.
                for amp in range(len(rpix_top)):
                    # Compute the ramp average for each pixel in this amplifier.
                    ra_top = np.mean(rpix_top[amp], axis=0)
                    ra_bot = np.mean(rpix_bot[amp], axis=0)
                    ramp_avg_top.append(ra_top)
                    ramp_avg_bot.append(ra_bot)

                # Loop over the vertical reference pixels to calculate their
                # ramp averages.
                for amp in range(len(rpix_side_fixed)):
                    # Compute the ramp average for each pixel in this amplifier.
                    ra_side = np.mean(rpix_side_fixed[amp], axis=0)
                    ramp_avg_side.append(ra_side)

                # Loop over the frames in the ramp.
                for frm in range(nframes):
                    # Create a list for the coefficients for this frame.
                    coeffs = []

                    # Loop over the amplifiers in this frame.
                    for amp in range(namps):
                        # Subtract the ramp average arrays for this amplifier.
                        rpix_top[amp][frm] -= ramp_avg_top[amp]
                        rpix_bot[amp][frm] -= ramp_avg_bot[amp]

                        # Compute the linear regression coefficients for this
                        # amplifier. Use the time array as the X axis and the
                        # ramp-averaged horizontal reference pixels as the
                        # Y axis.
                        rpix_horizontal = np.vstack([rpix_top[amp][frm], rpix_bot[amp][frm]]).flatten()
                        time_horizontal = np.vstack([time_top, time_bot]).flatten()

                        pfit = np.polyfit(time_horizontal, rpix_horizontal, deg=1)
                        coeffs.append(pfit)

                    # Loop over the amplifiers in this frame with the vertical
                    # reference pixels.
                    for amp in range(len(rpix_side_fixed)):
                        # Subtract the ramp average array for this amplifier.
                        rpix_side_fixed[amp][frm] -= ramp_avg_side[amp]

                    # Row average the top, side, and bottom reference pixels in
                    # the first amplifier. Stack them together into one array.
                    row_avg_top_left = np.mean(rpix_top[0][frm], axis=1, keepdims=True)
                    row_avg_bot_left = np.mean(rpix_bot[0][frm], axis=1, keepdims=True)
                    row_avg_side_left = np.mean(rpix_side_fixed[0][frm], axis=1, keepdims=True)
                    row_avg_left = np.vstack([row_avg_top_left, row_avg_side_left, row_avg_bot_left])

                    # Repeat for the last amplifier.
                    row_avg_top_right = np.mean(rpix_top[-1][frm], axis=1, keepdims=True)
                    row_avg_bot_right = np.mean(rpix_bot[-1][frm], axis=1, keepdims=True)
                    row_avg_side_right = np.mean(rpix_side_fixed[-1][frm], axis=1, keepdims=True)
                    row_avg_right = np.vstack([row_avg_top_right, row_avg_side_right, row_avg_bot_right])

                    # Repeat for the time array.
                    time_row_avg_top = np.mean(time_top, axis=1, keepdims=True)
                    time_row_avg_bot = np.mean(time_bot, axis=1, keepdims=True)
                    time_row_avg_side = np.mean(time_side, axis=1, keepdims=True)
                    time_row_avg = np.vstack([time_row_avg_top, time_row_avg_side, time_row_avg_bot])

                    # Construct a high frequency function for the first and
                    # last amplifiers. Average the two functions, row-wise.
                    highfreq_left = row_avg_left - coeffs[0][1] - (coeffs[0][0]*time_row_avg)
                    highfreq_right = row_avg_right - coeffs[-1][1] - (coeffs[-1][0]*time_row_avg)
                    highfreq = (highfreq_left + highfreq_right) / 2
                    ny, nx = highfreq.shape

                    # Now begin the FFT smoothing.
                    delt = 10E-6 * amp_shape[1]
                    smoothed_data = self._smoothFFT(highfreq, delt).reshape((ny, nx))

                    # Make the list for holding the amplifier correctors. Then
                    # loop over the amplifiers and create the corrector function
                    # for each amplifier. Each corrector should be shape
                    # (2048, 1). Replicate this to be (2048, 64) (the size of
                    # one amp). Then assemble the corrector array.
                    correctors = []
                    for amp in range(namps):
                        corr = coeffs[amp][1] + (coeffs[amp][0] * time_row_avg) + smoothed_data
                        corr = np.repeat(corr, amp_shape[1], 1)
                        correctors.append(corr)
                    corrector_array = np.hstack(correctors)

                    # Restructure and trim the corrector array to the size of
                    # the data array.
                    # This only satisfies full FOV imaging mode.
                    top_left = corrector_array[4:512, :4]       # vertical refpixels
                    bot_left = corrector_array[1536:2044, :4]   # vertical refpixels
                    top_right = corrector_array[4:512, -4:]     # vertical refpixels
                    bot_right = corrector_array[1536:2044, -4:] # vertical refpixels

                    top_zero = np.zeros((4,4)) # zero padding top vertical refpixels
                    mid_zero = np.zeros((8,4)) # zero padding middle vertical refpixels
                    bot_zero = np.zeros((4,4)) # zero padding bottom vertical refpixels

                    left_ref = np.concatenate([top_zero, top_left, mid_zero, bot_left, bot_zero], axis=0)
                    right_ref = np.concatenate([top_zero, top_right, mid_zero, bot_right, bot_zero], axis=0)

                    top_ref = corrector_array[:4, :]    # horizontal refpixels
                    bot_ref = corrector_array[-4:, :]   # horizontal refpixels

                    data = corrector_array[512:1536, :] # active pixels and middle vertical refpixels

                    corrector_array = np.concatenate([top_ref, data, bot_ref], axis=0)
                    corrector_array = np.concatenate([left_ref, corrector_array, right_ref], axis=1)

                    assert ext.data[frm].shape == corrector_array.shape

                    # Subtract the corrector frame from the data in the
                    # extension.
                    ext.data[frm] -= corrector_array

            # Update the filename.
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
            refpixt_kw = ad._keyword_for('ref_sec_top')
            refpixb_kw = ad._keyword_for('ref_sec_bot')
            refpixs_kw = ad._keyword_for('ref_sec_side')

            all_datasecs = ad.data_section()

            for ext, datasec in zip(ad, all_datasecs):
                # Trim SCI, VAR, DQ arrays
                #ext.reset(ext.nddata[:, datasec.y1:datasec.y2, datasec.x1:datasec.x2])

                # And OBJMASK (if it exists)
                #if hasattr(ext, 'OBJMASK'):
                #    ext.OBJMASK = ext.OBJMASK[:, datasec.y1:datasec.y2, datasec.x1:datasec.x2]

                # Temporarily use this to trim the arrays
                ext.data = ext.data[:, datasec.y1:datasec.y2, datasec.x1:datasec.x2]
                ext.variance = ext.variance[:, datasec.y1:datasec.y2, datasec.x1:datasec.x2]
                ext.mask = ext.mask[:, datasec.y1:datasec.y2, datasec.x1:datasec.x2]

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
                for amp in range(1,100):
                    if f'{darksec_kw}{amp}' in ext.hdr:
                        del ext.hdr[f'{darksec_kw}{amp}']
                    if f'{refsect_kw}{amp}' in ext.hdr:
                        del ext.hdr[f'{refsect_kw}{amp}']
                    if f'{refsecb_kw}{amp}' in ext.hdr:
                        del ext.hdr[f'{refsecb_kw}{amp}']
                    if f'{refsecs_kw}{amp}' in ext.hdr:
                        del ext.hdr[f'{refsecs_kw}{amp}']

            # Update the filename.
            ad.update_filename(suffix=sfx, strip=True)
        
        return adinputs

    def _smoothFFT(self, data, delt, first_deriv=False, second_deriv=False):
        """
        Optimal smoothing algorithm

        Smoothing algorithm to perform optimal filtering of the vertical
        reference pixels to reduce 1/f noise (horizontal stripes), based on the
        Kosarev and Pantos algorithm. This assumes that the data to be 
        filtered/smoothed has been sampled evenly.

        This smoothing algorithm was copied from pyNRC, which was adapted from
        M. Robberto's IDL code, optimal_smooth_fft.pro. No changes have been
        made to the pyNRC code.

        This code was accessed between October and December 2020. 

        pyNRC:
            https://github.com/JarronL/pynrc
        This code specifically:
            https://github.com/JarronL/pynrc/blob/master/pynrc/reduce/ref_pixels.py
        M. Robberto's code:
            http://www.stsci.edu/~robberto/Main/Software/IDL4pipeline/

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
            Fltr_Spectrum = np.zeros_like(Dat_P,dtype=np.complex)
            # make the filter complex
            i1 = 1; n2 = int((N-1)/2)+1; i2 = i1+n2
            FltrCoef = LOPT[i1:].astype(np.complex)
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