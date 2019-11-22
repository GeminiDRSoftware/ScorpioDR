"""
Recipes available to data with tags ['SCORPIO', 'CAL', 'ARC', 'NIR'].

Default recipe is set to "makeProcessedArc".
"""
recipe_tags = set(['SCORPIO', 'NIR', 'SPECT', 'LS', 'ARC'])

def makeProcessedArc(p):
    """
    This recipe performs the standardization and corrections needed to convert
    the raw input arc images into a reduced arc, with a wavelength solution
    attached to it.  This output processed arc is stored on disk using
    storeProcessedArc and has a name equal to the name of the first input arc
    image with "_arc.fits" appended.

    Parameters
    ----------
    p : Primitives object
        A primitive set matching the recipe_tags.
    """

    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.darkCorrect()  # TBD
    p.flatCorrect()  # TBD
    p.determineWavelengthSolution()
    p.determineDistortion()
    p.storeProcessedArc()
    p.writeOutputs()
    return

_default = makeProcessedArc


