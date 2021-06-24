#
#                                                                       DRAGONS
#
#                                            primitives_scorpio_image_nearIR.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt

from geminidr.core import Image, Photometry
from .primitives_scorpio_nearIR import ScorpioNearIR
from . import parameters_scorpio_nearIR_image

from recipe_system.utils.decorators import parameter_override

import numpy as np
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioNearIRImage(ScorpioNearIR, Image, Photometry):
    """
    This class contains primitives that applies to all Scorpio near-IR
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'IMAGE', 'NIR'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioNearIRImage, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_nearIR_image)

    def referencePixelCorrect(self, adinputs=None, **params):
        # Need to do checks for image size, like Full, Windowed, etc.
        # Need license info
        # Add citation to pyNRC and Massimo Robberto
        # Add further authoring info (me)
        # Code goes here

        """
        This primitive is used to correct a SCORPIO NIR image's noise using its reference pixels. SCORPIO's reference pixels are laid out very specifically as even the "full" image size is smaller than the detector size.

        The code in in this function

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """

        for ad in adinputs:
            # Get the list of keyword values for all reference pixel sections. 
            ref_kwds = ad.refpix_section()

            for ext in ad:
                # Convert any saturated or bad pixels to NaNs so they are not
                # used in any calculations.
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.saturated))] == np.nan
                ext.data[np.where(np.bitwise_and(ext.mask, DQ.bad_pixel))] == np.nan

                # Extract the horizontal reference pixel keyword values for 
                # this extension.
                rpix_top_kw = ref_kwds[ext]['top']
                rpix_bot_kw = ref_kwds[ext]['bottom']

                # Get the horizontal reference pixel data sections from the 
                # data plane. rpix_top and rpix_bot are lists of 3D arrays
                # (amplifiers) with shape (nframes, nrows, ncols).
                rpix_top = [ext.data[sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_top_kw]
                rpix_bot = [ext.data[sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_bot_kw]

                # Exctract the vertical reference pixel keyword values for
                # this extension.
                rpix_side_kw = ref_kwds[ext]['side']

                # Get the vertical reference pixel data sections from the 
                # data plane. rpix_side is a list of 3D arrays with shape
                # (nframes, nrows, ncols).
                rpix_side = [ext.data[sec.y1:sec.y2, sec.x1:sec.x2] for sec in rpix_side_kw]

                # Restructure the vertical reference pixel data sections into 
                # 2 sections (left and right). rpix_side_fixed is a list of 
                # 3D arrays (amplifiers) with shape (nframes, nrows, ncols).
                # As of June 2021, this is only verified for full FOV mode.
                rpix_left, rpix_right = [], []
                for sec in range(len(rpix_side)):
                    if sec % 2 == 0:
                        rpix_right.append(rpix_side[sec])
                    else:
                        rpix_left.append(rpix_side[sec])
                rpix_left = np.hstack(rpix_left)
                rpix_right = np.hstack(rpix_right)
                rpix_side_fixed = [rpix_side_left, rpix_side_right]







    def _smoothFFT():
        """
        Optimal smoothing algorithm

        Smoothing algorithm to perform optimal filtering of the vertical
        reference pixels to reduce 1/f noise (horizontal stripes), based on the
        Kosarev and Pantos algorithm. This assumes that the data to be 
        filtered/smoothed has been sampled evenly.

        This smooting algorithm was copied from pyNRC, which is based on
        M. Robberto's IDL code. No changes were made to the code.

        This code was accessed between October and December 2020. 

        pyNRC:
            https://github.com/JarronL/pynrc
        M. Robberto's code:
            http://www.stsci.edu/~robberto/Main/Software/IDL4pipeline/

        Parameters
        ----------


        Returns
        -------

        """



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