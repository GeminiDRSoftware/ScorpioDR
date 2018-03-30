__all__ = ['AstroDataOctocam']

from astrodata import factory
from gemini_instruments import addInstrumentFilterWavelengths
from .adclass import AstroDataOctocam
from .lookup import filter_wavelengths

factory.addClass(AstroDataOctocam)

addInstrumentFilterWavelengths('OCTOCAM', filter_wavelengths)

