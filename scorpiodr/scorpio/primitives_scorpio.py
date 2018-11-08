#
#                                                                       DRAGONS
#
#                                                         primitives_octocam.py
# ------------------------------------------------------------------------------

from geminidr.gemini.primitives_gemini import Gemini

from . import parameters_octocam

from .lookups import timestamp_keywords as octocam_stamps

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class OCTOCAM(Gemini):
    """
    This class inherits from the level above.  Any primitives specific
    to OCTOCAM and that applies to all channels can go here.  Primitives
    set specific to visible or nearIR, imaging or spectroscopy will go
    into subclasses inheriting OCTOCAM.
    """

    tagset = set()  # Cannot be assigned as a class

    def __init__(self, adinputs, **kwargs):
        super(OCTOCAM, self).__init__(adinputs, **kwargs)
        self._param_update(parameters_octocam)
        # Add OCTOCAM-specific timestamp keywords
        self.timestamp_keys.update(octocam_stamps.timestamp_keys)
