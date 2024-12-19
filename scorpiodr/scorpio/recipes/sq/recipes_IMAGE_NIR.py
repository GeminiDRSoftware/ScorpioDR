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
    p.addVAR(read_noise=True, poisson_noise=True)
    p.nonlinearityCorrect()
    p.referencePixelsCorrect()
    p.flagCosmicRaysFromNDRs()
    p.calculateSignalByRegression()

    p.darkCorrect()
    p.flatCorrect()
    p.stackIntegrations()

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
    p.adjustWCSToReference()
    p.resampleToCommonFrame()
    p.stackFrames()
    return

_default = reduce


def makeSkyFlat(p):
    """
    This recipe makes a flatfield image from a series of dithered sky images.

    Parameters
    ----------
    p : PrimitivesBASE object
        A primitive set matching the recipe_tags.
    """
    p.prepare()
    p.addDQ()
    p.ADUToElectrons()
    p.addVAR(read_noise=True, poisson_noise=True)
    p.nonlinearityCorrect()
    p.darkCorrect()
    # Make a "fastsky" by combining frames
    p.stackFrames(operation='median', scale=True, outstream='fastsky')
    p.normalizeFlat(stream='fastsky')
    p.thresholdFlatfield(stream='fastsky')
    # Flatfield with the fastsky and find objects
    p.flatCorrect(flat=p.streams['fastsky'][0], outstream='flattened')
    p.detectSources(stream='flattened')
    p.dilateObjectMask(dilation=10, stream='flattened')
    p.addObjectMaskToDQ(stream='flattened')
    p.transferAttribute(source='flattened', attribute='mask')
    p.stackFrames(operation='mean', scale=True, reject_method="minmax", nlow=0, nhigh=1)
    p.normalizeFlat()
    p.thresholdFlatfield()
    p.storeProcessedFlat(force=True)
    return
