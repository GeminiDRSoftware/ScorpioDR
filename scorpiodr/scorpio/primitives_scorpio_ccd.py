#
#                                                                       DRAGONS
#
#                                               primitives_scorpio_image_ccd.py
# ------------------------------------------------------------------------------

from geminidr.gemini.lookups import DQ_definitions as DQ
from gempy.gemini import gemini_tools as gt

from geminidr.core import CCD
from .primitives_scorpio import Scorpio
from . import parameters_scorpio_ccd

from recipe_system.utils.decorators import parameter_override
# ------------------------------------------------------------------------------

@parameter_override
class ScorpioCCD(Scorpio, CCD):
    """
    This class contains primitives that applies to all Scorpio optical
    imaging data.
    """

    tagset = set(['GEMINI', 'SCORPIO', 'CCD'])

    def __init__(self, adinputs, **kwargs):
        super(ScorpioCCD, self).__init__(adinputs, **kwargs)
        self.inst_lookups = 'scorpiodr.scorpio.lookups'
        self._param_update(parameters_scorpio_ccd)

    def trimOverscan(self, adinputs=None, suffix=None):
        """
        The trimOverscan primitive trims the overscan region from the input
        AstroData object and updates the headers.

        Parameters
        ----------
        suffix: str
            suffix to be added to output files
        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        for ad in adinputs:
            if ad.phu.get(timestamp_key) is not None:
                log.warning('No changes will be made to {}, since it has '
                            'already been processed by trimOverscan'.
                            format(ad.filename))
                continue

            ad = gt.trim_to_data_section(ad,
                                    keyword_comments=self.keyword_comments)
            # Tidy up additional header keywords that are too complicated
            # for gt.trim_to_data_section() to handle
            del ad.hdr['OVRSECS*']
            del ad.hdr['OVRSECP*']
            kw = ad._keyword_for('data_section')
            del ad.hdr[f'{kw}?*']

            # Set keyword, timestamp, and update filename
            ad.phu.set('TRIMMED', 'yes', self.keyword_comments['TRIMMED'])
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)
            ad.update_filename(suffix=suffix, strip=True)
        return adinputs

    def myNewPrimitive(self, adinputs=None, **params):
        """
        Description...

        Parameters
        ----------
        suffix: str
            suffix to be added to output files

        param2: blah
            blah, blah

        Returns
        -------

        """
        log = self.log
        log.debug(gt.log_message("primitive", self.myself(), "starting"))
        timestamp_key = self.timestamp_keys[self.myself()]

        # Get params out
        param2 = params['param2']

        # Initialize the list of output AstroData objects
        # It is also possible to modify adinputs in place.
        adoutputs = []

        for ad in adinputs:

            # Do whatever checks on the input are necessary, for example:
            # Check whether this primitive as been run already.
            if ad.phu.get(timestamp_key):
                log.warning("No changes will be made to {}, since it has"
                            "already been processed by myNewPrimitive".
                            format(ad.filename))
                continue

            # -----------------------
            # DR algorithm goes here
            # -----------------------

            # Timestamp
            gt.mark_history(ad, primname=self.myself(), keyword=timestamp_key)

            adoutputs.append(ad_out)

        return adoutputs