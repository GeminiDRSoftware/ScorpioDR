#                                                                       DRAGONS
#
#                                                  primitives_scorpio_bundle.py
# ------------------------------------------------------------------------------
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_bundle

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioBundle(Scorpio):
    """
    This is the class contains primitives that apply to Scorpio data that
    still contain all the channels.
    """
    tagset = set(["GEMINI", "SCORPIO", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioBundle, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_scorpio_bundle)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be verified.  I (KL) think that because the
        # data is all bundled up, the number of extension can be anything.
        # The default is 1 (in primitive_standardize).  So this function
        # here should override that.   Factor of 4 or 8, might be all we
        # can check.
        return len(ad) in [1, 4, 8]

    # The channel splitting primitive should go here
    # def separateChannels()