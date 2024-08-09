# This parameter file contains the parameters related to the primitives
# defined in the primitives_scorpio.py file.

from gempy.library import config

class stackIntegrationsConfig(config.Config):
    suffix = config.Field("Filename suffix", str, "_integrationsStacked", optional=True)


