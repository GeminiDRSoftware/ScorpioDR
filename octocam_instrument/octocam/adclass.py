
from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.common import Sextion
from gemini_instruments.gemini import AstroDataGemini

class AstroDataOctocam(AstroDataGemini):

    # single keyword mapping. add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict()

    @staticmethod
    def _matches_data(source):
        return source[0].header.get('INSTRUME', '').upper() == 'OCTOCAM'

    # ---------------
    # Tag definitions
    #----------------
    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['OCTOCAM'])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('OBSTYPE') == 'DARK':
            return TagSet(['DARK'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_arc(self):
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['ARC', 'CAL'])

    @astro_data_tag
    def _tag_bias(self):
        if self.phu.get('OBSTYPE') == 'BIAS':
            return TagSet(['BIAS', 'CAL'], blocks=['IMAGE', 'SPECT'])

    # More tags needs to be added by the OCTOCAM DR team
    # At this time, Gemini DR expects the following tags to be implemented.
    #    IMAGING, LS (for longslit), BUNDLE, FLAT, TWILIGHT, GCALFLAT.
    #    All type of flats must also be CAL and FLAT.
    #    CCD, NIR
    # Also probably needed:
    #    NODANDSHUFFLE, HIFREQ (High time resolution)

    # ----------------------
    # Descriptor definitions
    # ----------------------

    @astro_data_descriptor
    def overscan_section(self, pretty=False):
        """
        Returns the overscan (or bias) section.  If pretty is False, a
        tuple of 0-based coordinates is returned with format (x1, x2, y1, y2).
        If pretty is True, a keyword value is returned without parsing as a
        string.  In this format, the coordinates are generally 1-based.
        One tuple or string is return per extension/array.  If more than one
        array, the tuples/strings are return in a list.  Otherwise, the
        section is returned as a tuple or a string.

        Parameters
        ----------
        pretty : bool
         If True, return the formatted string found in the header.

        Returns
        -------
        tuple of integers or list of tuples
            Position of the overscan section using Python slice values.

        string or list of strings
            Position of the overscan section using an IRAF section
            format (1-based).
        """
        return self._parse_section('BIASSEC', pretty)

    # Obviously if BIASSEC is not the keyword used for OCTOCAM change that
    # in the example above.

    # For a list of the expected descriptors see the appendix in the
    # Astrodata User Manual.


