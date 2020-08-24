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

    __keyword_dict = dict(channel='CHANNEL',
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

    # Obviously if BIASSEC is not the keyword used for Scorpio change that
    # in the example above.

    # For a list of the expected descriptors see the appendix in the
    # Astrodata User Manual.


    # Modified 2020-08-21
    @astro_data_descriptor
    def amp_total_area(self, pretty=False):
        """
        Returns the total array size, including data and overscan for each amplifier in each extension.

        Returns
        -------
        list/(tuple or string)
            the total amp areas
        """

        # Get the channel and number of amps
        channel, namps = self._get_channel_amp_count()

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
        Returns a list of list of names of the amplifiers in the extensions
        or a single list of names if called on a single-extension slice.

        Returns
        -------
        list/list
            names of the amplifiers of the arrays
        """

        channel, namps = self._get_channel_amp_count()

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

    @astro_data_descriptor
    def channel(self):
        """
        Returns the channel name.

        Returns
        -------
        string
            Channel color band.
        """

        return self.phu.get(self._keyword_for('channel'), 1)

        '''
        if self.is_single:
            return self.phu.get(self._keyword_for('channel'), 1)
        else:
            return [self.phu.get(self._keyword_for('channel'), 1)]
        '''

    @astro_data_descriptor
    def gain(self):
        """
        Returns the gain (electrons/ADU) for each amplifier in each extension.
        
        Returns
        -------
        list
            Gains used for the observation
        """

        channel, namps = self._get_channel_amp_count()
        
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

    @astro_data_descriptor
    def read_noise(self):
        """
        Returns the read noise (electrons) for each amplifier in each extension.

        Returns
        -------
        list
            Read noise present in this observeration.
        """

        channel, namps = self._get_channel_amp_count()
        
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

    def _get_channel_amp_count(self):
        channel = self.channel().upper()
        namps = None

        if channel in ['G','R','I','Z']:
            namps = 4
        if channel in ['Y','J','H','K']:
            namps = 32

        return channel, namps


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