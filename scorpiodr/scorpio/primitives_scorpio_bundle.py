#                                                                       DRAGONS
#
#                                                  primitives_octocam_bundle.py
# ------------------------------------------------------------------------------
from .primitives_octocam import OCTOCAM
from . import parameters_octocam_bundle

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class OCTOCAMBundle(OCTOCAM):
    """
    This is the class contains primitives that apply to OCTOCAM data that
    still contain all the channels.
    """
    tagset = set(["GEMINI", "OCTOCAM", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(OCTOCAMBundle, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_octocam_bundle)

    # The channel splitting primitive should go here
    # def separateChannels()