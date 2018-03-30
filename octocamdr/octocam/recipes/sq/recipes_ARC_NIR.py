"""
Recipes available to data with tags ['OCTOCAM', 'CAL', 'ARC', 'NIR'].

Default recipe is set to "makeProcessedArc".
"""
recipe_tags = set(['GHOST', 'CAL', 'ARC', 'NIR'])

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
    p.addVAR(read_noise=True)
    p.nonlinearityCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    # whatever is needed
    #p.darkCorrect()
    #p.extractSpectrum()
    #p.fitWavelength()
    p.storeProcessedArc()
    return

default = makeProcessedDark