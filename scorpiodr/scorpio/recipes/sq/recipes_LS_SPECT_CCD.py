"""
Recipes available to data with tags ['SCORPIO', 'SPECT', 'CCD'].
Default is "reduce".
"""
recipe_tags = set(['SCORPIO', 'SPECT', 'CCD'])

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
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    #p.darkCorrect()
    #p.scatteredLightCorrect()
    #p.applyWavelengthSolution()    # depends on Gemini's algorithm.
    #p.QECorrect()                  # depends on Gemini's algorithm.
    p.flatCorrect()
    p.rejectCosmicRays()   # TDB
    p.distortionCorrect()
    p.findSourceApertures()
    p.skyCorrectFromSlit()
    p.resampleToCommonFrame()
    p.stackFrames()
    p.findSourceApertures()
    p.traceApertures()
    p.extract1DSpectra()
    p.fluxCalibrate()
    p.linearizeSpectra()   # TBD
    p.writeOutputs()
    return

_default = reduce


def reduceStandard(p):
    """
    todo: add docstring

    Parameters
    ----------
    p : :class:`geminidr.core.primitives_gmos_longslit.GMOSLongslit`

    """
    p.prepare()
    p.addDQ(static_bpm=None)
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.flatCorrect()
    p.distortionCorrect()
    p.findSourceApertures(max_apertures=1)
    p.skyCorrectFromSlit()
    p.resampleToCommonFrame()
    p.stackFrames()
    p.traceApertures()
    p.extract1DSpectra()
    p.linearizeSpectra()   # TDB if needed.
    p.calculateSensitivity()
    p.writeOutputs()