#                                                                       DRAGONS
#
#                                                  primitives_scorpio_bundle.py
# ------------------------------------------------------------------------------
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_bundle

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioBundle(Scorpio):
    """
    This is the class contains primitives that apply to Scorpio data that
    still contain all the channels.
    """
    tagset = set(["GEMINI", "SCORPIO", "BUNDLE"])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioBundle, self).__init__(adinputs, **kwargs)
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

    # The channel splitting primitive should go here
    # def separateChannels()

    # Added by Landon Gelman, Jun 29, 2020
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
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has "
                            "already been processed by {}".format(
                            ad.filename, self.myself()))
                continue

            log.stdinfo("Unbundling {}:".format(ad.filename))

            # We just want all extensions written to their own file, so do all.
            # TODO : decide what the suffix should be for new files.
            #        This should be added to _write_newfile()
            for ext in ad:
                _write_newfile(ext, ad, log)


def _write_newfile(extns, base, log):
    """
    Helper function to write sub-files out from a MEF bundle.

    Parameters
    ----------
    extns : iterable of :any:`astrodata.Astrodata`
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
        If the ``extns`` parameter is :any:`None`, or empty
    """

    assert extens and len(extens) > 0

    # Copy the PHU from the astrodata base object
    nf = astrodata.create(copy.deepcopy(base.phu))

    # Make changes to the PHU
        # BUNDLE to "F"
        # NEXTEND to 1

    try:
        nf.phu.set('NEXTEND', 1)
        nf.phu.set('BUNDLE', 'F')
    except KeyError:
        pass
    
    # Append the extension to the newly created ad object.
    for x in extns:
        if (x.hdr.get('NAXIS') > 0) and (x.data.size > 0):
            newfile.append(x)
    nf.filename = base.filename
    nf.update_filename(


# TODO
# All of GHOST's _get_hdr_values() and _get_common_hdr_values() can be skipped. They're using it as a check to make sure all the extensions contain the same value for the same keyword. But in Scorpio's case, any common keyword will already be in the primary header.
# Copy the primary header to each individual extension
# Get the extension and append to new astrodata object (along with original PHU)
# Get the suffix to be added on to the filename
# Write the new FITS file
# Repeat with the next extension from the original MEF