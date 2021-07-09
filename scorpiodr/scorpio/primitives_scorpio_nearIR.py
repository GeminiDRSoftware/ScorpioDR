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
        Two-Point Difference method for finding outliers in a 3-D ramp data array.
        The scheme used in this variation of the method uses numpy array methods
        to compute first-differences and find the max outlier in each pixel while
        still working in the full 3-D data array. This makes detection of the first
        outlier very fast. We then iterate pixel-by-pixel over only those pixels
        that are already known to contain an outlier, to look for any additional
        outliers and set the appropriate DQ mask for all outliers in the pixel.
        This is MUCH faster than doing all the work on a pixel-by-pixel basis.

        This function has been modified from its original location, here:
        https://github.com/spacetelescope/stcal
        Accessed February 2021-.

        The method used in this function is based on the method by Anderson & Gordon, 2011, which can be found here: https://iopscience.iop.org/article/10.1086/662593


        """

    def referencePixelCorrect(self, adinputs=None, **params):
        """
        This primitive is used to correct a SCORPIO NIR image's noise using its
        reference pixels. SCORPIO's reference pixels are laid out very
        specifically as even the "full" image size is smaller than the detector
        size.

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

    def _smoothFFT(self, data, delt, first_deriv=False, second_deriv=False):
        """
        Optimal smoothing algorithm

        Smoothing algorithm to perform optimal filtering of the vertical
        reference pixels to reduce 1/f noise (horizontal stripes), based on the
        Kosarev and Pantos algorithm. This assumes that the data to be 
        filtered/smoothed has been sampled evenly.

        This smooting algorithm was copied from pyNRC, which was adapted from
        M. Robberto's IDL code. No changes have been made to the pyNRC code.

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