#
#                                                                       DRAGONS
#
#                                                  primitives_scorpio_bundle.py
# ------------------------------------------------------------------------------
import copy

import astrodata
from gempy.gemini import gemini_tools as gt
from recipe_system.utils.decorators import parameter_override

from .primitives_scorpio import Scorpio
from . import parameters_scorpio_bundle
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioBundle(Scorpio):
    """
    This class contains primitives that apply to Scorpio data that
    still contain all the channels.
    """
    tagset = set(["GEMINI", "SCORPIO", "BUNDLE"])

    def _initialize(self, adinputs, **kwargs):
        #super(ScorpioBundle, self).__init__(adinputs, **kwargs)
        super()._initialize(adinputs, **kwargs)
        self._param_update(parameters_scorpio_bundle)

    @staticmethod
    def _has_valid_extensions(ad):
        """ Check that the AD has a valid number of extensions. """

        # this needs to be verified.  I (KL) think that because the
        # data is all bundled up, the number of extension can be anything.
        # The default is 1 (in primitive_standardize).  So this function
        # here should override that.   Factor of 4 or 8, might be all we
        # can check.
        return len(ad) in [1, 4, 8]

    def separateChannels(self, adinputs=None, **params):
        """
        Break a Scorpio observation bundle into individual exposures.

        This primitive separates a Scorpo observation bundle into individual 
        exposures. The number of files created will equal the number of 
        extensions in the original bundle.
        """

        # Create string that will be the suffix appended to the filename for de-bundled images.
        # GHOST uses the individual CAMERA and EXPID keywords.

        log = self.log
        log.debug(gt.log_message('primitive', self.myself(), 'starting'))
        #timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            """
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by {}".format(
                            ad.filename, self.myself()))
                continue
            """

            log.stdinfo("Unbundling {}:".format(ad.filename))

            # We just want all extensions written to their own file, so do all.
            # TODO : decide what the suffix should be for new files.
            #        This should be added to _write_newfile()
            for idx, ext in enumerate(ad):
                _write_newfile(idx, ext, ad, log)

        return []


def _write_newfile(enum, ext, base, log):
    """
    Helper function to write sub-files out from a MEF bundle.

    Parameters
    ----------
    ext : iterable of :any:`astrodata.Astrodata`
        AstroData extensions to be appended to the new file
    base : :any:`astrodata.Astrodata`
        Original AstroData instance to base the new file on. The file's primary 
        header unit will be used, as will the base of the filename.
    log : AstroData logging object
        Log recording actions. This should be the log in use in the calling 
        primitive.

    Raises
    ------
    AssertionError
        If the ``ext`` parameter is :any:`None`, or empty
    """

    assert ext and len(ext) > 0

    # Copy the PHU from the astrodata base object
    nf = astrodata.create(copy.deepcopy(base.phu))

    # Make changes to the PHU
        # BUNDLE to "F"
        # NEXTEND to 1
        # Add CHANNEL keyword from extension HDR to new file PHU

    try:
        nf.phu.set('NEXTEND', 1)
        nf.phu.set('BUNDLE', 'F')
        nf.phu.set('CHANNEL', ext.hdr['CHANNEL'])
        nf.phu.set('XOFFSET', ext.hdr['XOFFSET'])
        nf.phu.set('YOFFSET', ext.hdr['YOFFSET'])
        nf.phu.set('POFFSET', ext.hdr['POFFSET'])
        nf.phu.set('QOFFSET', ext.hdr['QOFFSET'])
        nf.phu.set('RAOFFSET', ext.hdr['RAOFFSET'])
        nf.phu.set('DECOFFSE', ext.hdr['DECOFFSE'])
        #nf.phu.set('HA', ext.hdr['HA'])
        #nf.phu.set('ELEVATIO', ext.hdr['ELEVATIO'])
        #nf.phu.set('AZIMUTH', ext.hdr['AZIMUTH'])
        #nf.phu.set('PAR_ANG', ext.hdr['PAR_AND'])
        #nf.phu.set('PA', ext.hdr['PA'])
    except KeyError:
        pass
    else:
    #    del ext.hdr['CHANNEL']
        del ext.hdr['XOFFSET']
        del ext.hdr['YOFFSET']
        del ext.hdr['POFFSET']
        del ext.hdr['QOFFSET']
        del ext.hdr['RAOFFSET']
        del ext.hdr['DECOFFSE']
        #del ext.hdr['HA']
        #del ext.hdr['ELEVATIO']
        #del ext.hdr['AZIMUTH']
        #del ext.hdr['PAR_ANG']
        #del ext.hdr['PA']

    
    # Append the extension to the newly created ad object.
    for x in ext:
        if (x.hdr.get('NAXIS') > 0) and (x.data.size > 0):
            nf.append(x)

    channel = ext.hdr['CHANNEL']
    nf.phu.set('FILE', ('_').join(base.filename.split('_')[:-1]) + '_{}_{}.fits'.format(enum+1, channel))
    nf.filename = ('_').join(base.filename.split('_')[:-1]) + '_{}_{}.fits'.format(enum+1, channel)

    log.stdinfo("    Writing {}".format(nf.filename))
    nf.write(overwrite=True)


# TODO: Fix the suffix to be added on to the filename
