__all__ = ['AstroDataScorpio']

from astrodata import factory
from gemini_instruments.gemini import addInstrumentFilterWavelengths
from .adclass import AstroDataScorpio
from .lookup import filter_wavelengths

factory.addClass(AstroDataScorpio)

addInstrumentFilterWavelengths('Scorpio', filter_wavelengths)

