"""
Recipes available to data with tags ['SCORPIO', 'IMAGE', 'CCD'].

Default is "reduce".
"""
recipe_tags = set(['SCORPIO', 'IMAGE', 'CCD'])

def reduce(p):
    """
    This recipe performs the standardization and corrections needed to
    convert the raw input science images into a stacked image.

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
    p.darkCorrect()   # TBD if truly needed
    p.flatCorrect()
    p.makeFringe()
    p.fringeCorrect()
    p.detectSources()
    p.alignAndStack()
    #p.detectSources()  # might be good to create the catalog from the stack.
    p.fluxCalibrate()  # magic to find calibrator.  adds a zeropoint.
    return

default = reduce


def makeProcessedFringe(p):
    """
    This recipe creates a fringe frame from the inputs files. The output
    is stored on disk using storeProcessedFringe and has a name equal
    to the name of the first input bias image with "_fringe.fits" appended.

    Parameters
    ----------
    p : PrimitivesBASEE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.addDQ()
    p.addVAR(read_noise=True)
    p.overscanCorrect()
    p.biasCorrect()
    p.ADUToElectrons()
    p.addVAR(poisson_noise=True)
    p.flatCorrect()
    p.makeFringeFrame()
    p.storeProcessedFringe()
    return