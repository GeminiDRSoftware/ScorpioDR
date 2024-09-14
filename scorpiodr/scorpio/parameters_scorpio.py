# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio.py file.

from geminidr.core import parameters_visualize
from gempy.library import config

class stackIntegrationsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_integrationsStacked", optional=True)


class displayConfig(parameters_visualize.displayConfig):
    remove_bias = config.Field("Remove estimated bias level before displaying?",
                               bool, True)
    integration = config.RangeField("Integration to display (1-indexed)", int,
                                    min=1, default=None, optional=True)

    def setDefaults(self):
        self.tile = False
