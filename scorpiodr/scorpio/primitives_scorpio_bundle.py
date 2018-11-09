#                                                                       DRAGONS
#
#                                                  primitives_scorpio_bundle.py
# ------------------------------------------------------------------------------
from .primitives_scorpio import Scorpio
from . import parameters_octocam_bundle

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
        self._param_update(parameters_octocam_bundle)

    # The channel splitting primitive should go here
    # def separateChannels()