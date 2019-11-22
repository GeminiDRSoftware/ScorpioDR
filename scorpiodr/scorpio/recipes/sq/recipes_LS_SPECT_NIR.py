"""
Recipes available to data with tags ['SCORPIO', 'SPECT', 'NIR'].
Default is "reduce".
"""
recipe_tags = set(['SCORPIO', 'SPECT', 'NIR'])

def reduce(p):
    """
    This recipe process optical spectrum up to and including alignment and
    stacking.  A single stacked extracted and calibrated spectrum is produced.


    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.nonlinearityCorrect()
    p.darkCorrect()
    #p.applyWavelengthSolution()    # depends on Gemini's algorithm.
    #p.QECorrect()                  # depends on Gemini's algorithm.
    p.flatCorrect()
    p.separateSky()
    p.associateSky()
    p.skyCorrect()
    p.distortionCorrect()
    p.findSourceApertures()
    p.traceApertures()
    p.extract1DSpectra()
    p.linearizeSpectra()
    p.telluricCorrect()
    p.stackSpectra()
    p.writeOutputs()
    return

_default = reduce


# are the tellurics identified as CAL?  If so, move to its own recipe file.

def reduceTelluric(p):
    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.nonlinearityCorrect()
    p.darkCorrect()
    # p.applyWavelengthSolution()    # depends on Gemini's algorithm.
    # p.QECorrect()                  # depends on Gemini's algorithm.
    p.flatCorrect()
    p.separateSky()
    p.associateSky()
    p.skyCorrect()
    p.distortionCorrect()
    p.findSourceApertures()
    p.traceApertures()
    p.extract1DSpectra()
    p.linearizeSpectra()
    p.stackSpectra()
    p.writeOutputs()
    return
