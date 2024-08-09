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