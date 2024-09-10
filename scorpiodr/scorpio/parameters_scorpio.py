# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio.py file.

from geminidr.core import parameters_visualize
from gempy.library import config

class stackIntegrationsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_integrationsStacked", optional=True)


class displayConfig(parameters_visualize.displayConfig):
    def setDefaults(self):
        self.tile = False
