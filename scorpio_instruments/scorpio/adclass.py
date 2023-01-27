from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.gemini import AstroDataGemini
from gemini_instruments.common import Section

import numpy as np

from . import lookup


def tuple_to_section(sec, pretty=False):
    return sec.asIRAFSection() if pretty else sec


class AstroDataScorpio(AstroDataGemini):

    # single keyword mapping. add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict(array_name='ARRNAM',
                          array_section='ARRSEC',
                          channel='CHANNEL',
                          dark_section='DARKSEC',
                          data_section='DATSEC',
                          detector_name='DETECTOR',
                          read_noise='RDNOIS',
                          ref_sec_top='REFSCT',
                          ref_sec_bot='REFSCB',
                          ref_sec_side='REFSCS',
                          )

    @staticmethod
    def _matches_data(source):
        return source[0].header.get('INSTRUME', '').upper() == 'SCORPIO'

    # ---------------
    # Tag definitions
    #----------------
    @astro_data_tag
    def _tag_instrument(self):
        return TagSet(['SCORPIO'])

    @astro_data_tag
    def _is_bundle(self):
        if self.phu.get('BUNDLE') == 'T':
            return TagSet(['BUNDLE'])
        else:
            return TagSet(blocks=['BUNDLE'])

    @astro_data_tag
    def _tag_dark(self):
        if self.phu.get('OBSTYPE') == 'DARK':
            return TagSet(['DARK', 'CAL'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_arc(self):
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['ARC', 'CAL'])

    @astro_data_tag
    def _tag_bias(self):
        if self.phu.get('OBSTYPE') == 'BIAS':
            return TagSet(['BIAS', 'CAL', 'CCD'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('OBSTYPE') == 'FLAT':
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_is_ccd(self):
        if (self.phu.get('BUNDLE') == 'F') and (self.phu.get('CHANNEL').upper() in ['G','R','I','Z']):
            return TagSet(['CCD'], blocks=['NIR'])

    @astro_data_tag
    def _tag_is_nir(self):
        if (self.phu.get('BUNDLE') == 'F') and (self.phu.get('CHANNEL').upper() in ['Y','J','H','K']):
            return TagSet(['NIR'], blocks=['CCD'])

    @astro_data_tag
    def _flat_type(self):
        bun = self.phu.get('BUNDLE')
        obj = self.phu.get('OBJECT', '').upper()
        shut = self.phu.get('GCALSHUT', '').upper()
        if bun == 'F' and obj == 'GCALFLAT' and shut == 'OPEN':
            return TagSet(['LAMPON', 'NIR'], blocks=['CCD'])
        if bun == 'F' and obj == 'GCALFLAT' and shut == 'CLOSED':
            return TagSet(['LAMPOFF', 'NIR'], blocks=['CCD'])

    # More tags needs to be added by the Scorpio DR team
    # At this time, Gemini DR expects the following tags to be implemented.
    #    IMAGING, LS (for longslit), BUNDLE, FLAT, TWILIGHT, GCALFLAT.
    #    All type of flats must also be CAL and FLAT.
    #    CCD, NIR
    # Also probably needed:
    #    NODANDSHUFFLE, HIFREQ (High time resolution)

    @astro_data_descriptor
    def array_name(self):
        """
        Returns the name for each amplifier array per extension. Because Scorpio
        has multiple amplifiers per extension, this returns a list of strings 
        per extension. If the method is called on a single slice, the names are 
        returned in a list of strings. Otherwise the names are returned in a 
        list of lists of strings.

        Returns
        -------
        list of str / list of list of str
            Names of the amplifiers of the arrays.
        """

        values = []
        keyword = self._keyword_for('array_name')
        for amp in range(1,100):
            value = self.hdr.get(f'{keyword}{amp}')
            if self.is_single:
                if value is None:
                    break
                values.append(value)
            else:
                if value[0] is None:
                    break
                values.append(value[0])
        if self.is_single:
            return values
        else:
            return [values]

    @astro_data_descriptor
    def array_section(self, pretty=False):
        """
        Returns the section covered by the array(s) relative to the detector 
        frame. For example, this can be the position of multiple amps read 
        within a CCD. If pretty is False, a tuple of 0-based coordinates 
        is returned with format (x1, x2, y1, y2). If pretty is True, a keyword 
        value is returned without parsing as a string. In this format, the 
        coordinates are generally 1-based.

        In the case of Scorpio, each extension returns either a list of tuples
        (one per amplifier), or a single string containing a 1-indexed section
        per amp, joined with commas. If the method is called on a single slice,
        the sections are returned as tuples or a single string.

        Parameters
        ----------
        pretty : bool
            If True, return a list of formatted strings found in the header.

        Returns
        -------
        list of tuple of integers or list of list of tuples
            Positions of arrays in extension(s) using Python slice values.

        str/list of str
            Position of arrays in extension(s) using a 1-based section format.
        """
        arrsec = self._build_section_lists(self._keyword_for('array_section'))
        if self.is_single:
            return (tuple_to_section(arrsec, pretty=pretty)
                    if isinstance(arrsec, Section) else
                    (",".join(tuple_to_section(sec, pretty=True) for sec in arrsec)
                     if pretty else arrsec))

        return [tuple_to_section(asec, pretty=pretty)
                if isinstance(asec, Section) else
                (",".join(tuple_to_section(sec, pretty=True) for sec in asec)
                 if pretty else asec) for asec in arrsec]

    @astro_data_descriptor
    def channel(self):
        """
        Returns the channel name. Returns a string if the Scorpio file is 
        de-bundled or a list if the Scorpio file is a bundle.

        Returns
        -------
        list of string/string
            Channel color band.
        """
        return self.phu.get(self._keyword_for('channel'))

    @astro_data_descriptor
    def data_section(self, pretty=False):
        """
        Returns the rectangular section that includes the pixels that would be
        exposed to light.  If pretty is False, a tuple of 0-based coordinates
        is returned with format (x1, x2, y1, y2).  If pretty is True, a keyword
        value is returned without parsing as a string.  In this format, the
        coordinates are generally 1-based.

        One tuple or string is return per extension/array, in a list. If the
        method is called on a single slice, the section is returned as a tuple
        or a string.

        For SCORPIO, this descriptor assumes that the individual amps sections
        are contiguous, and therefore a single section is returned per extension

        Parameters
        ----------
        pretty : bool
            if True, return a 1-indexed string representation

        Returns
        -------
        tuple/str or list of tuple/str
            location of the pixels exposed to light
        """
        datasec = self._build_section_lists(self._keyword_for('data_section'))
        if self.is_single:
            if isinstance(datasec, list):
                datasec = Section(x1=min(s.x1 for s in datasec), x2=max(s.x2 for s in datasec),
                                  y1=min(s.y1 for s in datasec), y2=max(s.y2 for s in datasec))
            return tuple_to_section(datasec, pretty=pretty)

        sections = [Section(x1=min(s.x1 for s in dsec), x2=max(s.x2 for s in dsec),
                            y1=min(s.y1 for s in dsec), y2=max(s.y2 for s in dsec))
                    if isinstance(dsec, list) else dsec for dsec in datasec]
        return [tuple_to_section(sec, pretty=pretty) for sec in sections]

    @astro_data_descriptor
    def gain(self):
        """
        Returns the gain (electrons/ADU) for each amplifier in each extension. 
        Because Scorpio has multiple amplifiers per extension, this returns a 
        list of floats per extension.
        
        Returns
        -------
        list of floats / list of list of floats
            Gains used for the observation.
        """

        values = []
        keyword = self._keyword_for('gain')
        for amp in range(1, 100):
            value = self.hdr.get(f'{keyword}{amp}')
            if self.is_single:
                if value is None:
                    break
                values.append(value)
            else:
                if value[0] is None:
                    break
                values.append(value[0])

        if self.is_single:
            return values
        else:
            return [values]

    @astro_data_descriptor
    def overscan_section(self, pretty=False):
        """
        Returns the overscan (or bias) sections.  If pretty is False, each
        section is returned as a tuple of 0-based coordinates with format
        (x1, x2, y1, y2). If pretty is True, a keyword value is returned
        without parsing as a string. In this format, the coordinates are
        generally 1-based. The descriptor for SCORPIO returns a dict keyed
        by 'serial' and 'parallel' with each value being either a single
        section in the format dictated by "pretty" (for a single extension)
        or a list of such sections, one per extension.

        Parameters
        ----------
        pretty : bool
         If True, return the formatted string found in the header.

        Returns
        -------
        dict
        """
        try:
            overscan_dict = {'serial': self._build_section_lists('OVRSECS', pretty=pretty)}
        except KeyError:
            # Something for IR arrays?
            return None if self.is_single else [None] * len(self)
        else:
            overscan_dict['parallel'] = self._build_section_lists('OVRSECP', pretty=pretty)
            return overscan_dict

    @astro_data_descriptor
    def read_noise(self):
        """
        Returns the read noise (electrons) for each amplifier in each extension. 
        Because Scorpio has multiple amplifiers per extension, this returns a 
        list of floats per extension.
        
        Returns
        -------
        list of floats / list of list of floats
            Read noised present in the observation.
        """

        values = []
        keyword = self._keyword_for('read_noise')
        for amp in range(1, 100):
            value = self.hdr.get(f'{keyword}{amp}')
            if self.is_single:
                if value is None:
                    break
                values.append(value)
            else:
                if value[0] is None:
                    break
                values.append(value[0])

        if self.is_single:
            return values
        else:
            return [values]

    @astro_data_descriptor
    def refpix_section(self, pretty=False):
        topsec = self._build_section_lists(self._keyword_for('ref_sec_top'))
        botsec = self._build_section_lists(self._keyword_for('ref_sec_bot'))
        sidesec = self._build_section_lists(self._keyword_for('ref_sec_side'))

        if self.is_single:
            top = (tuple_to_section(topsec, pretty=pretty) 
                   if isinstance(topsec, Section) else
                   (",".join(tuple_to_section(sec, pretty=True) for sec in topsec)
                    if pretty else topsec))
            bot = (tuple_to_section(botsec, pretty=pretty) 
                   if isinstance(botsec, Section) else
                   (",".join(tuple_to_section(sec, pretty=True) for sec in botsec)
                    if pretty else botsec))
            side = (tuple_to_section(sidesec, pretty=pretty) 
                    if isinstance(sidesec, Section) else
                    (",".join(tuple_to_section(sec, pretty=True) for sec in sidesec)
                     if pretty else sidesec))
            return ({'top':top, 'bottom':bot, 'side':side})

        top = [tuple_to_section(tsec, pretty=pretty)
               if isinstance(tsec, Section) else
               (",".join(tuple_to_section(sec, pretty=True) for sec in tsec)
                if pretty else tsec) for tsec in topsec]
        bot = [tuple_to_section(bsec, pretty=pretty)
               if isinstance(bsec, Section) else
               (",".join(tuple_to_section(sec, pretty=True) for sec in bsec)
                if pretty else bsec) for bsec in botsec]
        side = [tuple_to_section(ssec, pretty=pretty)
               if isinstance(ssec, Section) else
               (",".join(tuple_to_section(sec, pretty=True) for sec in ssec)
                if pretty else ssec) for ssec in sidesec]
        return ({'top':top, 'bottom':bot, 'side':side})

    def _build_section_lists(self, keyword, pretty=False):
        # See if there is only one keyword without a number
        sec = self._parse_section(keyword, pretty=pretty)
        if not (sec is None or not self.is_single and sec.count(None) == len(sec)):
            return sec

        # OK, find and resort the keywords
        sections = []
        for amp in range(1, 100):
            sec = self._parse_section(f'{keyword}{amp}', pretty=pretty)
            if sec is None or not self.is_single and sec.count(None) == len(sec):
                break
            sections.append(sec)
        if amp == 1:
            raise KeyError(f"Keywords {keyword} and {keyword}1 not found")
        if self.is_single:
            return sections
        return [[sec[i] for sec in sections if sec[i] is not None]
                for i in range(len(self))]

