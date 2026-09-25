import math

from astrodata import astro_data_tag, astro_data_descriptor, returns_list, TagSet
from gemini_instruments import gmu
from gemini_instruments.gemini import AstroDataGemini, get_specphot_name
from gemini_instruments.common import Section, build_group_id

from . import lookup

def tuple_to_section(sec, pretty=False):
    return sec.asIRAFsection() if pretty else sec


class AstroDataScorpio(AstroDataGemini):

    # single keyword mapping. add only the ones that are different
    # from what's already defined in AstroDataGemini.

    __keyword_dict = dict(array_name='ARRNAM',
                          array_section='ARRSEC',
                          dark_section='UNILSEC',
                          data_section='DATSEC',
                          detector_name='DETECTOR',
                          disperser='GRATING',
                          filter='FILTER',
                          read_noise='RDNOIS',
                          slit='SLITSIZE',
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

    def _tag_image_or_spect(self):
        mode = self.phu.get(self._keyword_for('observation_mode'), '').upper()
        if mode.startswith('IM'):
            return 'IMAGE'
        elif mode.startswith('SPEC') or mode == 'N&S':
            return 'SPECT'

    def _tag_is_ccd(self):
        return self.phu.get('CHANNEL', '') in ['g','r','i','z']

    def _tag_is_nir(self):
        return self.phu.get('CHANNEL', '') in ['Y','J','H','Ks']

    def _tag_is_bias(self):
        return self.phu.get('OBSTYPE') == 'BIAS'

    def _tag_is_bpm(self):
        return self.phu.get('OBSTYPE') == 'BPM' or 'BPMASK' in self.phu

    def _tag_is_dark(self):
        return self.phu.get('OBSTYPE') == 'DARK'

    @astro_data_tag
    def _tag_arc(self):
        if self.phu.get('OBSTYPE') == 'ARC':
            return TagSet(['ARC', 'CAL'])

    @astro_data_tag
    def _tag_bias(self):
        if self._tag_is_bias():
            return TagSet(['BIAS', 'CAL', 'CCD'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_dark(self):
        if self._tag_is_dark():
            return TagSet(['DARK', 'CAL'], blocks=['IMAGE', 'SPECT'])

    @astro_data_tag
    def _tag_flat(self):
        if self.phu.get('OBSTYPE') == 'FLAT':
            return TagSet(['FLAT', 'CAL'])

    @astro_data_tag
    def _tag_standard(self):
        if (
            self._tag_image_or_spect() == 'SPECT' and
            self.phu.get('OBSTYPE') == 'OBJECT' and
            self.phu.get('OBSCLASS') in ('partnerCal', 'nightCal') and
            (self._tag_is_nir() or get_specphot_name(self))
        ):
            return TagSet(['STANDARD', 'CAL'])

    @astro_data_tag
    def _tag_twilight(self):
        if self.phu.get('OBJECT', '').upper() == 'TWILIGHT':
            # Twilight flats are of OBSTYPE == OBJECT, meaning that the generic
            # FLAT tag won't be triggered. Add it explicitly
            return TagSet(
                ['TWILIGHT', 'CAL',
                 'SLITILLUM' if self._tag_image_or_spect() == 'SPECT' else
                 'FLAT']
            )

    @astro_data_tag
    def _tag_ccd(self):
        if self._tag_is_ccd():
            return TagSet(['CCD'], blocks=['NIR'])

    @astro_data_tag
    def _tag_nir(self):
        if self._tag_is_nir():
            return TagSet(['NIR'], blocks=['CCD'])

    @astro_data_tag
    def _type_mode(self):  # overrides the like-named gemini tag method
        modes = {'IMAGE', 'SPECT'}
        mode = self._tag_image_or_spect()
        if mode in modes:
            other = list(modes - {mode})
            return TagSet([mode], blocks=other)

    @astro_data_tag
    def _tag_ls(self):
        if self._tag_is_bias() or self._tag_is_dark() or self._tag_is_bpm():
            return

        if str(self.phu.get('SLITSIZE', '')).endswith('arcsec'):
            return TagSet(['LS'])

    @astro_data_tag
    def _tag_nodandshift(self):
        if self.phu.get('NODMODE', '').lower() == 'yes':
            return TagSet(['NODANDSHIFT'])

    @astro_data_tag
    def _flat_type(self):
        obj = self.phu.get('OBJECT', '').upper()
        shut = self.phu.get('GCALSHUT', '').upper()
        if obj == 'GCALFLAT' and shut == 'OPEN':
            return TagSet(['LAMPON']) #, 'NIR'], blocks=['CCD'])
        if obj == 'GCALFLAT' and shut == 'CLOSED':
            return TagSet(['LAMPOFF']) #, 'NIR'], blocks=['CCD'])

    # More tags needs to be added by the Scorpio DR team
    # At this time, Gemini DR expects the following tags to be implemented.
    #    IMAGING, LS (for longslit), BUNDLE, FLAT, TWILIGHT, GCALFLAT.
    #    All type of flats must also be CAL and FLAT.
    #    CCD, NIR
    # Also probably needed:
    #    NODANDSHUFFLE, HIFREQ (High time resolution)

    @astro_data_descriptor
    def amp_read_area(self):
        """
        Returns a list of amplifier read areas, one per amp, made by combining
        the amplifier name and detector section (nested within a length-1 list
        of extensions when called on the parent AstroData instance).

        Returns
        -------
        list[Section|str] | list[list[Section|str]]
            read_area of each extension
        """
        ampname = self.array_name()
        arrsec = self.array_section(pretty=True)
        # Combine the amp name(s) and detector section(s)
        if self.is_single:
            return ["'{}':{}".format(a, s) if a and s else None
                    for a, s in zip(ampname, arrsec)]
        else:
            return [["'{}':{}".format(a, s) if a and s else None
                     for a, s in zip(amps, secs)]
                    for amps, secs in zip(ampname, arrsec)]

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

        In the case of Scorpio, each extension returns a list of either Section
        objects or (for pretty=True) strings containing 1-indexed sections, one
        per amplifier. When this method is called on a top-level AstroData
        instance, those are nested in an outer (usually length-1) list of
        extensions.

        Parameters
        ----------
        pretty : bool
            If True, return a list of formatted strings found in the header.

        Returns
        -------
        list of tuple of integers or list of list of tuples
            Positions of arrays in extension(s) using Python slice values.

        list[str] | list[list[str]]
            Position of arrays in extension(s) using a 1-based section format.
        """
        arrsec = self._build_section_lists(self._keyword_for('array_section'))
        if self.is_single:
            return (tuple_to_section(arrsec, pretty=pretty)
                    if isinstance(arrsec, Section) else
                    (list(tuple_to_section(sec, pretty=True) for sec in arrsec)
                     if pretty else arrsec))

        return [tuple_to_section(asec, pretty=pretty)
                if isinstance(asec, Section) else
                (list(tuple_to_section(sec, pretty=True) for sec in asec)
                 if pretty else asec) for asec in arrsec]

    @astro_data_descriptor
    def camera(self, stripID=False, pretty=False):
        """
        Returns the name of the camera.

        The 'stripID' & 'pretty' options currently have no effect, as SCORPIO
        doesn't append a component ID anyway.

        Returns
        -------
        str
            The name of the camera (eg. 'g').

        """
        return self.channel()

    @astro_data_descriptor
    @gmu.return_requested_units(input_units="um")
    def central_wavelength(self):
        """
        Returns the central wavelength for a spectrum (in m by default)

        Returns
        -------
        float
            The central wavelength setting
        """

        val = self.phu.get(self._keyword_for('central_wavelength'), None)

        if val is None:
            chan = self.channel() or self.filter_name(pretty=True)
            val = lookup.central_wavelengths.get(chan)

        return float(val) if val else None

    @astro_data_descriptor
    def channel(self):
        """
        Returns the channel name.

        Returns
        -------
        string
            Channel color band.
        """
        return self.phu.get('CHANNEL')

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
    def detector_roi_setting(self):
        """
        Returns the ROI setting.

        Returns
        -------
        str
            Name of the ROI setting used: "Full Frame" or "Window"
            (or "Undefined" if not recognized).
        """
        roi_dict = lookup.ROI_settings
        roi_settings = set()
        for ext in self:
            roi = ext.detector_section()
            roi_setting = 'Undefined'
            for s in roi_dict:
                roi_tuple = (roi.y1, roi.y2, roi.x1, roi.x2)
                if roi_tuple in roi_dict[s]:
                    roi_setting = s
            roi_settings.add(roi_setting)
        return roi_settings.pop() if len(roi_settings)==1 else 'Undefined'

    @astro_data_descriptor
    def detector_x_bin(self):
        """
        Returns the detector binning in the x-direction

        Returns
        -------
        int
            The detector binning
        """

        def _get_xbin(b):
            try:
                return int(b.split()[0])
            except (AttributeError, ValueError):
                return None

        binning = self.hdr.get('CCDSUM')
        if self.is_single:
            return _get_xbin(binning) if 'CCD' in self.tags else 1
        else:
            if 'CCD' in self.tags:
                xbin_list = [_get_xbin(b) for b in binning]
            else:
                xbin_list = [1 for ext in self]
            # Check list is single-valued
            return xbin_list[0] if xbin_list[1:] == xbin_list[:-1] else None

    @astro_data_descriptor
    def detector_y_bin(self):
        """
        Returns the detector binning in the y-direction

        Returns
        -------
        int
            The detector binning
        """

        def _get_ybin(b):
            try:
                return int(b.split()[1])
            except (AttributeError, ValueError, IndexError):
                return None

        binning = self.hdr.get('CCDSUM')
        if self.is_single:
            return _get_ybin(binning) if 'CCD' in self.tags else 1
        else:
            if 'CCD' in self.tags:
                ybin_list = [_get_ybin(b) for b in binning]
            else:
                ybin_list = [1 for ext in self]
            # Check list is single-valued
            return ybin_list[0] if ybin_list[1:] == ybin_list[:-1] else None

    @astro_data_descriptor
    @gmu.return_requested_units()
    def dispersion(self):
        """
        Returns the dispersion in nm per pixel as a list (one value per
        extension) or a float if used on a single-extension slice. It is
        possible to control the units of wavelength using the input arguments.

        Returns
        -------
        list/float
            The dispersion(s) in m/pixel
        """

        chan = self.channel()
        dispersion = lookup.dispersions.get(chan)
        xbin = self.detector_x_bin()

        if dispersion and xbin and 'SPECT' in self.tags:
            dispersion *= xbin
        else:
            dispersion = None

        if not self.is_single:
            dispersion = [dispersion] * len(self)

        return dispersion

    @returns_list
    @astro_data_descriptor
    def dispersion_axis(self):
        """
        Returns the axis along which the light is dispersed.

        Returns
        -------
        (list of) int (1)
            Dispersion axis.
       """
        return 1

    @astro_data_descriptor
    def filter_name(self, stripID=False, pretty=False, keepID=False):
        """
        Returns the name of the filter(s) used.  The component ID can be
        removed with either 'stripID' or 'pretty'.

        Parameters
        ----------
        stripID : bool
            If True, removes the component ID and returns only the name of
            the filter.
        pretty : bool
            Same as for stripID.  Pretty here does not do anything more.

        Returns
        -------
        str
            The name of the filter with or without the component ID.
        """
        return self._may_remove_component(self._keyword_for('filter'),
                                          stripID, pretty)

    @astro_data_descriptor
    def focal_plane_mask(self, stripID=False, pretty=False):
        """
        Returns the name of the focal plane mask.

        Parameters
        ----------
        stripID : bool
            Doesn't actually do anything.
        pretty : bool
            Same as for stripID

        Returns
        -------
        str
            The name of the focal plane mask
        """
        # Does this need changing to incorporate both rails, like GNIRS,
        # and/or strip component numbers?
        mask = self.slit()
        return 'Imaging' if mask == 'None' else mask  # same convention as GMOS

    @astro_data_descriptor
    def gain(self):
        """
        Returns the gain (electrons/ADU) for each amplifier in each extension. 
        Because Scorpio has multiple amplifiers per extension, this returns a 
        list of floats per extension for raw data.
        
        Returns
        -------
        list of floats / list of list of floats
            Gains used for the observation.
        """

        values = []
        keyword = self._keyword_for('gain')

        for ext in self:
            extval = []
            try:
                extval = ext.hdr[keyword]  # once GAIN is set, overrides GAINn
            except KeyError:
                for amp in range(1, 100):
                    value = ext.hdr.get(f'{keyword}{amp}')
                    if value is None:
                        break
                    extval.append(value)
            values.append(extval)

        if self.is_single:
            return values[0]
        else:
            return values

    @astro_data_descriptor
    def gain_setting(self):
        """
        Returns the gain setting of the observation (which is fixed for
        SCORPIO).

        Returns
        -------
        str
            Gain setting
        """
        return 'Standard'

    @astro_data_descriptor
    def group_id(self):
        """
        Returns a string representing a group of data that are compatible
        with each other. This is used when stacking, for example.

        SCORPIO uses the detector binning, amp_read_area and
        read_speed_setting. Data other than biases & darks have the pretty
        version of the filter name included, while darks have the exposure
        time. Spectroscopic data have the grism and slit. Nod-and-shift cals
        have the charge shuffling parameters. Science data & standards include
        the observation_id.

        Returns
        -------
        str
            A group ID for compatible data.
        """
        tags = self.tags

        # Things needed for all observations
        id_descriptor_list_all = ['detector_x_bin', 'detector_y_bin',
                                  'read_mode', 'amp_read_area']

        if 'BIAS' in tags:
            id_descriptor_list = []
        elif 'DARK' in tags:
            id_descriptor_list = ['exposure_time']
        else:
            if 'STANDARD' in tags or 'CAL' not in tags:
                id_descriptor_list = ['observation_id', 'filter_name']
            else:
                id_descriptor_list = ['filter_name']

            # Matching disperser as well as filter is redundant for SCORPIO,
            # but would be needed if matching components without "pretty":
            if 'SPECT' in tags:
                id_descriptor_list.extend(['disperser', 'focal_plane_mask'])

        # Match charge shuffling parameters for cals, since they get stacked
        # early on and this will affect the pixel flat & charge smearing.
        # Science data usually get stacked later on and we can just match their
        # OBSID, allowing stacking of irregular nod sequences if that's what
        # the PI indended (note that GMOS doesn't account for N&S at all here).
        if 'NODANDSHIFT' in tags and 'CAL' in tags and not 'STANDARD' in tags:
            id_descriptor_list.extend(['shuffle_pixels', 'nod_count'])

        # Add in all of the common descriptors required
        id_descriptor_list.extend(id_descriptor_list_all)

        return build_group_id(self, id_descriptor_list,
                              prettify=['filter_name', 'disperser'],  # why?
                              additional=None)

    @astro_data_descriptor
    def nod_count(self):
        """
        Returns a tuple with the number of integrations made in each
        of the nod-and-shift positions

        Returns
        -------
        tuple of int
            number of integrations in the A and B positions
        """
        try:
            return (int(self.phu['ANODCNT']), int(self.phu['BNODCNT']))
        except KeyError:
            return None

    @astro_data_descriptor
    def nod_offsets(self):
        """
        Returns a tuple with the offsets from the default telescope position
        of the A and B nod-and-shift positions (in arcseconds)

        Returns
        -------
        tuple of float
            offsets in arcseconds
        """
        try:
            ayoff = float(self.phu['NODAYOFF'])
            byoff = float(self.phu['NODBYOFF'])
            inport = self.phu['INPORT']
        except KeyError:
            return None

        if inport == 1:
            ayoff, byoff = -ayoff, -byoff  # TBC for SCORPIO IAA/detectors

        return (ayoff, byoff)

    @astro_data_descriptor
    def nominal_photometric_zeropoint(self):
        """
        Returns the nominal zeropoints (i.e., the magnitude corresponding to
        a pixel value of 1) for the extensions in an AD object.
        Zeropoints in table are for electrons, so subtract 2.5*log10(avg_gain)
        if the data are in ADU.

        Returns
        -------
        float/list
            zeropoint values, one per SCI extension
        """

        def _zpt(ext):
            det = ext.detector_name()
            try:
                zpt = lookup.nominal_zeropoints[
                    (det, ext.filter_name())  # allow exact component match
                ]
            except KeyError:
                zpt = lookup.nominal_zeropoints.get(
                    (det, ext.filter_name(pretty=True))
                )

            in_adu = ext.is_in_adu()
            gain = ext.gain()

            try:
                # Just use average gain across the quadrants if we still have
                # ADU; we'll get a slightly more accurate value at the final
                # pass with fully-reduced data:
                return zpt - (
                    2.5 * math.log10(math.fsum(gain)/len(gain))
                    if in_adu else 0
                )
            except (TypeError, ZeroDivisionError):
                return None

        if self.is_single:
            return _zpt(self)
        else:
            return [_zpt(ext) for ext in self]

    @astro_data_descriptor
    def non_linear_level(self):
        # temporary value for testing dragons compatibility before real data
        return self.saturation_level()

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
            try:
                overscan_dict['parallel'] = self._build_section_lists('OVRSECP', pretty=pretty)
            except KeyError:
                pass
            return overscan_dict

    @astro_data_descriptor
    def read_mode(self):
        """
        Returns a string describing the readout mode, which sets the
        readout speed

        Returns
        -------
        str
            read mode used
        """
        read_mode = self.phu.get('RDMODE')
        if read_mode == 0:
            return 'Standard'
        elif read_mode == 1:
            return 'Fast'
        else:
            return None

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
    def read_speed_setting(self):
        """
        Returns the setting for the readout speed (Standard or Fast)

        Returns
        -------
        str
            the setting for the readout speed
        """
        # The gain setting is fixed for SCORPIO and we derive the read speed
        # from the read mode because that's how it's defined in the headers.
        return self.read_mode()

    @astro_data_descriptor
    def refpix_section(self, pretty=False):
        useSidePixels = True if self.phu.get('OBSMODE') == "spect" or self.phu.get('OBSMODE') == 'image' and self.phu.get('IMTYPE') == "full" else False

        topsec = self._build_section_lists('REFSCT', pretty=pretty)
        botsec = self._build_section_lists('REFSCB', pretty=pretty)
        if useSidePixels:
            sidesec = self._build_section_lists('REFSCS', pretty=pretty)

        if self.is_single:
            top = (tuple_to_section(topsec, pretty=pretty) 
                   if isinstance(topsec, Section) else
                   (",".join(tuple_to_section(sec, pretty=True) for sec in topsec)
                    if pretty else topsec))
            bot = (tuple_to_section(botsec, pretty=pretty) 
                   if isinstance(botsec, Section) else
                   (",".join(tuple_to_section(sec, pretty=True) for sec in botsec)
                    if pretty else botsec))
            if useSidePixels:
                side = (tuple_to_section(sidesec, pretty=pretty) 
                        if isinstance(sidesec, Section) else
                        (",".join(tuple_to_section(sec, pretty=True) for sec in sidesec)
                         if pretty else sidesec))
                return ({'top':top, 'bottom':bot, 'side':side})
            return ({'top':top, 'bottom':bot})

        top = [tuple_to_section(tsec, pretty=pretty)
               if isinstance(tsec, Section) else
               (",".join(tuple_to_section(sec, pretty=True) for sec in tsec)
                if pretty else tsec) for tsec in topsec]
        bot = [tuple_to_section(bsec, pretty=pretty)
               if isinstance(bsec, Section) else
               (",".join(tuple_to_section(sec, pretty=True) for sec in bsec)
                if pretty else bsec) for bsec in botsec]
        if useSidePixels:
            side = [tuple_to_section(ssec, pretty=pretty)
                   if isinstance(ssec, Section) else
                   (",".join(tuple_to_section(sec, pretty=True) for sec in ssec)
                    if pretty else ssec) for ssec in sidesec]
            return ({'top':top, 'bottom':bot, 'side':side})
        return ({'top':top, 'bottom':bot})

    @astro_data_descriptor
    def saturation_level(self):
        # temporary value for testing dragons compatibility before real data
        level= 65535
        return level if self.is_single else [level for ext in self]

    @astro_data_descriptor
    def shuffle_pixels(self):
        """
        Returns the number of rows that the charge has been shuffled, in
        nod-and-shift data

        Returns
        -------
        int
            The number of rows by which the charge is shuffled
        """
        nodpix = self.phu.get('NODPIX')
        if 'NODANDSHIFT' in self.tags and nodpix is not None:
            return int(nodpix)  # Draft GDS config defines a string...

    @astro_data_descriptor
    def slit_width(self):
        """
        Returns the width of the slit in arcseconds

        Returns
        -------
        float/None
            the slit width in arcseconds
        """
        fpmask = self.slit()
        if fpmask and 'arcsec' in fpmask:
            return float(fpmask.replace('arcsec', ''))
        return None

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

