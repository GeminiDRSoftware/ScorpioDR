"""
Recipes available to data with tags ['SCORPIO', 'SPECT', 'CCD'].
Default is "reduceScience".
"""
recipe_tags = set(['SCORPIO', 'SPECT', 'CCD'])

def reduceScience(p):
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
    p.attachWavelengthSolution()
    #p.QECorrect()                  # depends on Gemini's algorithm.
    p.flatCorrect()
    p.flagCosmicRays()
    #p.applyDQPlane(replace_flags=8, inner=1.0, outer=5.0)
    p.distortionCorrect()
    p.findApertures()
    p.skyCorrectFromSlit()
    p.adjustWCSToReference()
    p.resampleToCommonFrame(conserve=True)  # default output_wave_scale="linear"
    p.scaleCountsToReference()
    p.stackFrames()
    p.findApertures()
    p.traceApertures()
    p.storeProcessedScience(suffix="_2D")
    p.extractSpectra()
    p.fluxCalibrate()
    p.storeProcessedScience(suffix="_1D")
    p.writeOutputs()
    return

_default = reduceScience


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
    p.attachWavelengthSolution()
    p.flatCorrect()
    p.distortionCorrect()
    p.findApertures(max_apertures=1)
    p.skyCorrectFromSlit()
    p.resampleToCommonFrame(conserve=True)  # default output_wave_scale="linear"
    p.scaleCountsToReference()
    p.stackFrames()
    p.traceApertures()
    p.extractSpectra()
    p.calculateSensitivity()
    p.storeProcessedStandard()
    p.writeOutputs()
