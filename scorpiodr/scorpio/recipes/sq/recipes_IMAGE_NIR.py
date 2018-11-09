"""
Recipes available to data with tags ['SCORPIO', 'IMAGE', 'NIR'].
Default is "reduce".
"""
recipe_tags = set(['SCORPIO', 'IMAGE', 'NIR'])

def reduce(p):
    """
    This recipe process near-IR image up to and including alignment and
    stacking.  A single stacked output image is produced.


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

    # Initial sky subtraction (pre-masking)
    p.separateSky()
    p.associateSky(stream='sky')
    p.skyCorrect(instream='sky', mask_objects=False, outstream='skysub')

    # mask sources in sky frames
    p.detectSources(stream='skysub')
    p.transferAttribute(stream='sky', source='skysub', attribute='OBJMASK')
    p.clearStream(stream='skysub')

    # proper sky subtraction with source masked.
    p.associateSky()
    p.skyCorrect(mask_objects=True)
    p.detectSources()
    p.alignAndStack()
    #p.detectSources()   # would be a good thing to have.
    p.fluxCalibrate()  # magic to find calibrator.  adds a zeropoint.
    return

default = reduce

# p.prepare()
# p.addDQ()
# p.ADUToElectrons()
# p.addVAR(read_noise=True, poisson_noise=True)
# p.nonlinearityCorrect()
# p.darkCorrect()
# p.flatCorrect()
# p.separateSky()
# p.associateSky()
# p.skyCorrect(mask_objects=False)
# p.alignAndStack()
# return
