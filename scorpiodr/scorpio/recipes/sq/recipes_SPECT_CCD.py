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
    p.darkCorrect()
    p.scatteredLightCorrect()
    #p.applyWavelengthSolution()
    p.QECorrect()
    p.flatCorrect()
    p.rejectCosmicRays()
    p.applyWavelengthSolution()  # transform
    p.skyCorrect()
    p.fluxCalibrate()
    p.stackSpectra()
    p.extractSpectrum()
    return

_default = reduce
