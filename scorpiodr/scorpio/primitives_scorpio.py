#
#                                                                       DRAGONS
#
#                                                         primitives_scorpio.py
# ------------------------------------------------------------------------------

from geminidr.gemini.primitives_gemini import Gemini

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
        self._param_update(parameters_scorpio)
        # Add Scorpio-specific timestamp keywords
        self.timestamp_keys.update(scorpio_stamps.timestamp_keys)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be verified. This would be for the debundled data
        # and inherited and used by the ScorpioImage and ScorpioSpect class.
        return len(ad) in [1]
