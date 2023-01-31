#
#                                                                       DRAGONS
#
#                                                         primitives_scorpio.py
# ------------------------------------------------------------------------------

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
            ad.update_filename(suffix, strip=True)
        return adinputs