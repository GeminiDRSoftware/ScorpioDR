"""
Recipes available to data with tags ['SCORPIO', 'SPECT', 'NIR'].
Default is "reduce".
"""
recipe_tags = set(['SCORPIO', 'SPECT', 'NIR'])

def reduce(p):
    """
    This recipe process near-IR spectrum up to and including alignment and
    stacking.  A single stacked extracted and calibrated spectrum is produced.


    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.nonlinearityCorrect()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.darkCorrect()
    p.flatCorrect()
    p.skyCorrect()  # ABBA etc, sky substraction
    p.stackSpectra()  # spatial shift of B to A, and stack.
    p.applyWavelengthSolution()  # transform, how about s-distortion?
    p.extractSpectrum()
    p.telluricCorrect()
    return

_default = reduce


# are the tellurics identified as CAL?  If so, move to its own recipe file.

def reduceTelluric(p):
    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.nonlinearityCorrect()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.darkCorrect()
    p.flatCorrect()
    p.skyCorrect()  # ABBA etc, sky substraction
    p.stackSpectra()  # spatial shift of B to A, and stack.
    p.applyWavelengthSolution()  # transform, how about s-distortion?
    p.extractSpectrum()
    return
