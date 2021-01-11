from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.common import Section
from gemini_instruments.gemini import AstroDataGemini
from gemini_instruments.common import section_to_tuple

import numpy as np

from . import lookup

class AstroDataScorpio(AstroDataGemini):

    # single keyword mapping. add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict(array_name='ARRNAM',
                          array_section='ARRSEC',
                          channel='CHANNEL',
                          data_section='DATSEC',
                          read_noise='RDNOIS',
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

    # More tags needs to be added by the Scorpio DR team
    # At this time, Gemini DR expects the following tags to be implemented.
    #    IMAGING, LS (for longslit), BUNDLE, FLAT, TWILIGHT, GCALFLAT.
    #    All type of flats must also be CAL and FLAT.
    #    CCD, NIR
    # Also probably needed:
    #    NODANDSHUFFLE, HIFREQ (High time resolution)

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('OBSTYPE') == 'FLAT':
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_is_ccd(self):
        #if ('BUNDLE' not in self.tags) and (self.phu.get('CHANNEL').upper() in ['G','R','I','Z']):
        if (self.phu.get('BUNDLE') == 'F') and (self.phu.get('CHANNEL').upper() in ['G','R','I','Z']):
            return TagSet(['CCD'], blocks=['NIR'])

    @astro_data_tag
    def _tag_is_nir(self):
        #if ('BUNDLE' not in self.tags) and (self.phu.get('CHANNEL').upper() in ['Y','J','H','K']):
        if (self.phu.get('BUNDLE') == 'F') and (self.phu.get('CHANNEL').upper() in ['Y','J','H','K']):
            return TagSet(['NIR'], blocks=['CCD'])

    # TODO: Check the work I did in ScorpioDR_old for GCAL tags

    @astro_data_tag
    def _flat_type(self):
        bun = self.phu.get('BUNDLE')
        obj = self.phu.get('OBJECT', '').upper()
        shut = self.phu.get('GCALSHUT', '').upper()
        if bun == 'F' and obj == 'GCALFLAT' and shut == 'OPEN':
            return TagSet(['LAMPON', 'NIR'], blocks=['CCD'])
        if bun == 'F' and obj == 'GCALFLAT' and shut == 'CLOSED':
            return TagSet(['LAMPOFF', 'NIR'], blocks=['CCD'])

    # Write our own _parse_section based on the method in the Gemini adclass
    # It should take into consideration the number of amplifiers in Scorpio
    """
    def _parse_section(self, keyword, pretty, namps):
        section_list = []
        try:
            value_filter = (str if pretty else section_to_tuple)
            process_fn = lambda x: (None if x is None else value_filter(x))
            section = []
            for _i in range(1,namps+1):
                section.append(self.hdr.get(keyword+'{}'.format(_i)))
            if self.is_single:
                return
    """

    # ----------------------
    # Descriptor definitions
    # ----------------------

    @astro_data_descriptor
    def detector_name(self):
        """
        Returns the name of the detector

        Returns
        -------
        str
            the detector name
        """
        return self.phu.get(self._keyword_for('detector_name'))


    # For a list of the expected descriptors see the appendix in the
    # Astrodata User Manual.


    # Modified 2020-08-27
    @astro_data_descriptor
    def amplifier_layout(self):
        """
        Returns a tuple of the multi-amplifier shape/layout when multiple 
        amplifiers exist on one extension. If called on more than one extension,
        the tuples are returned in a list. Otherwise the shape/layout is 
        returned as a tuple.

        Returns
        -------
        tuple of int / list of tuples of ints
            Shape or layout of multi-amplifier configuration per extension.
        """
        channel = self.channel().upper()
        namps = lookup.amplifier_count.get(channel)
        
        arrsec = self.array_section(pretty=False)

        if self.is_single:
            v_amps = 0
            h_amps = 0
            v_max = 0
            h_max = 0
            for arr in arrsec:
                if arr[1] > h_max:
                    h_amps += 1
                    h_max = arr[1]
                if arr[3] > v_max:
                    v_amps += 1
                    v_max = arr[3]
            amp_shape = (v_amps,h_amps)
            return amp_shape
        else:
            amp_shape = []
            for sec in arrsec:
                v_amps = 0
                h_amps = 0
                v_max = 0
                h_max = 0
                for arr in sec:
                    if arr[1] > h_max:
                        h_amps += 1
                        h_max = arr[1]
                    if arr[3] > v_max:
                        v_amps += 1
                        v_max = arr[3]
                amp_shape.append((v_amps,h_amps))
            return amp_shape

    @astro_data_descriptor
    def amp_total_area(self, pretty=False):
        """
        Returns the total array size, including data and overscan for each 
        amplifier in each extension.

        Returns
        -------
        list/(tuple or string)
            the total amp areas
        """

        # Get the channel and number of amps
        #channel, namps = self._get_channel_amp_count()
        channel = self.channel().upper()
        namps = lookup.amplifier_count.get(channel)

        # Create a list for the sections
        read_arrays = []

        # Loop over the namps to get all keywords
        for amp in range(1,namps+1):
            # Get values for keywords DATSECn, OVRSECSn, and OVRSECPn
            datsec = self._parse_section('DATSEC{}'.format(amp), pretty=False)
            ovrsecs = self._parse_section('OVRSECS{}'.format(amp), pretty=False)
            ovrsecp = self._parse_section('OVRSECP{}'.format(amp), pretty=False)

            # Concatenate values to new strings like [x1:x2,y1:y2]
            if self.is_single:
                x1 = np.min([datsec[0], ovrsecs[0], ovrsecp[0]])
                x2 = np.max([datsec[1], ovrsecs[1], ovrsecp[1]])
                y1 = np.min([datsec[2], ovrsecs[2], ovrsecp[2]])
                y2 = np.max([datsec[3], ovrsecs[3], ovrsecp[3]])
            else:
                x1 = np.min([datsec[0][0], ovrsecs[0][0], ovrsecp[0][0]])
                x2 = np.max([datsec[0][1], ovrsecs[0][1], ovrsecp[0][1]])
                y1 = np.min([datsec[0][2], ovrsecs[0][2], ovrsecp[0][2]])
                y2 = np.max([datsec[0][3], ovrsecs[0][3], ovrsecp[0][3]])

            amparr = "[{0}:{1},{2}:{3}]".format(x1+1, x2, y1+1, y2)
            
            value_filter = (str if pretty else section_to_tuple)
            process_fn = lambda x: (None if x is None else value_filter(x))

            read_arrays.append(process_fn(amparr))

        if self.is_single:
            return read_arrays
        else:
            return [read_arrays]

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

        channel = self.channel().upper()
        namps = lookup.amplifier_count.get(channel)
        keyword = self._keyword_for('array_name')

        name_list = []
        for amp in range(1,namps+1):
            name = self.hdr.get(keyword+'{}'.format(amp))
            if self.is_single:
                name_list.append(name)
            else:
                name_list.append(name[0])
        if self.is_single:
            return name_list
        else:
            return [name_list]

        '''
        names = []

        for amp in range(1,namps+1):
            name = self.hdr.get('ARRNAM{}'.format(amp))
            names.append(name)

        if self.is_single:
            return names
        else:
            return [names]
        '''
        '''
        names = []
        if self.is_single:
            for amp in range(1,namps+1):
                name = self.hdr.get('ARRNAM{}'.format(amp))
                names.append(name)
        else:
            for ext in range(len(self)):
                extnames = []
                for amp in range(1,namps+1):
                    name = self[ext].hdr.get('ARRNAM{}'.format(amp))
                    extnames.append(name)
                names.append(extnames)
        return names
        '''

    @astro_data_descriptor
    def array_section(self, pretty=False):
        """
        Returns the section covered by the array(s) relative to the detector 
        frame. For example, this can be the position of multiple amps read 
        within a CCD. If pretty is False, a tuple of 0-based coordinates 
        is returned with format (x1, x2, y1, y2). If pretty is True, a keyword 
        value is returned without parsing as a string. In this format, the 
        coordinates are generally 1-based.

        In the case of Scorpio, one list of tuples or strings is returned per 
        extension/array in a list. If the method is called on a single slice, 
        the sections are returned as tuples or strings in a list.

        Parameters
        ----------
        pretty : bool
            If True, return a list of formatted strings found in the header.

        Returns
        -------
        list of tuple of integers or list of list of tuples
            Positions of arrays in extension(s) using Python slice values.

        list of str/list of list of str
            Position of arrays in extension(s) using a 1-based section format.
        """
        arrsec = self._build_section_lists(self._keyword_for('array_section'))
        if self.is_single:
            section = Section(x1=min(s.x1 for s in arrsec), x2=max(s.x2 for s in arrsec),
                              y1=min(s.y1 for s in arrsec), y2=max(s.y2 for s in arrsec))
            if pretty:
                return f'[{section.x1+1}:{section.x2}:{section.y1+1}:{section.y2}]'
            return section

        section = [Section(x1=min(s.x1 for s in asec), x2=max(s.x2 for s in asec),
                           y1=min(s.y1 for s in asec), y2=max(s.y2 for s in asec))
                   for asec in arrsec]
        if pretty:
            return [f'[{s.x1+1}:{s.x2}:{s.y1+1}:{s.y2}]' for s in section]
        return section

    @astro_data_descriptor
    def detector_section(self, pretty=False):
        arrsec = self._build_section_lists(self._keyword_for('detector_section'))
        if self.is_single:
            section = Section(x1=min(s.x1 for s in arrsec), x2=max(s.x2 for s in arrsec),
                              y1=min(s.y1 for s in arrsec), y2=max(s.y2 for s in arrsec))
            if pretty:
                return f'[{section.x1+1}:{section.x2}:{section.y1+1}:{section.y2}]'
            return section

        section = [Section(x1=min(s.x1 for s in asec), x2=max(s.x2 for s in asec),
                           y1=min(s.y1 for s in asec), y2=max(s.y2 for s in asec))
                   for asec in arrsec]
        if pretty:
            return [f'[{s.x1+1}:{s.x2}:{s.y1+1}:{s.y2}]' for s in section]
        return section

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
        keyword = self._keyword_for('data_section')
        return self._build_section_lists(keyword, pretty=pretty)

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

        channel = self.channel().upper()
        namps = lookup.amplifier_count.get(channel)
        keyword = self._keyword_for('gain')

        gain_list = []
        for amp in range(1,namps+1):
            gain = self.hdr.get(keyword+'{}'.format(amp))
            if self.is_single:
                gain_list.append(gain)
            else:
                gain_list.append(gain[0])
        if self.is_single:
            return gain_list
        else:
            return [gain_list]

        """
        gain_list = []
        if self.is_single:
            for i in range(namps):
                try:
                    gain = self.hdr.get('GAIN{}'.format(i+1))
                except Exception:
                    break
                else:
                    gain_list.append(gain)
        else:
            for i in range(len(self)):
                ext_gain_list = []
                for j in range(namps):
                    try:
                        gain = self[i].hdr.get('GAIN{}'.format(j+1))
                    except Exception:
                        break
                    else:
                        ext_gain_list.append(gain)
                gain_list.append(ext_gain_list)
        
        return gain_list
        """

    @astro_data_descriptor
    def map_amplifier_to_layout(self, keyword):
        """
        Returns the amplifier-specific values (keyword based) in the 
        shape of the amplifier order (described in lookup.py) on the CCD 
        detector.

        Returns
        -------

        """
        channel = self.channel().upper()
        layout = self.amplifier_layout()
        order = lookup.amplifier_order.get(channel)
        
        amparr = []
        for amp in order:
            kw = self.hdr.get(keyword+'{}'.format(amp))
            if self.is_single:
                amparr.append(kw)
            else:
                amparr.append(kw[0])
        
        amparr = np.array([amparr])

        if self.is_single:
            amparr = amparr.reshape(layout)
            return amparr
        else:
            amparr = amparr.reshape(layout[0])
            return [amparr]

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

        channel = self.channel().upper()
        namps = lookup.amplifier_count.get(channel)
        keyword = self._keyword_for('read_noise')

        readnoise_list = []
        for amp in range(1,namps+1):
            readnoise = self.hdr.get(keyword+'{}'.format(amp))
            if self.is_single:
                readnoise_list.append(readnoise)
            else:
                readnoise_list.append(readnoise[0])
        if self.is_single:
            return readnoise_list
        else:
            return [readnoise_list]
        
        """
        rn_list = []
        if self.is_single:
            for i in range(namps):
                try:
                    rn = self.hdr.get('RDNOIS{}'.format(i+1))
                except Exception:
                    break
                else:
                    rn_list.append(rn)
        else:
            for i in range(len(self)):
                ext_rn_list = []
                for j in range(namps):
                    try:
                        rn = self[i].hdr.get('RDNOIS{}'.format(j+1))
                    except Exception:
                        break
                    else:
                        ext_rn_list.append(rn)
                rn_list.append(ext_rn_list)
        
        return rn_list
        """

    '''
    def _get_channel_amp_count(self):
        channel = self.channel().upper()
        namps = None

        if channel in ['G','R','I','Z']:
            namps = 4
        if channel in ['Y','J','H','K']:
            namps = 32

        return channel, namps
    '''
    '''
    @astro_data_descriptor
    def gain(self):
        """
        Return the averaged gain for the color channel.

        Returns
        -------
        float
            gain
        """
        channel = self.channel().upper()
        gain = lookup.array_properties.get('gain{}'.format(channel))
        return [gain]
    '''
    '''
    @astro_data_descriptor
    def read_noise(self):
        """
        Return the averaged read noise for the color channel.
        
        Returns
        -------
        float
            read noise
        """
        channel = self.channel().upper()
        readnoise = lookup.array_properties.get('read_noise{}'.format(channel))
        return [readnoise]
    '''

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
