"""
Recipes available to data with tags ['Scorpio', 'IMAGE', 'NIR'].
Default is "reduce".
"""
recipe_tags = set(['Scorpio', 'IMAGE', 'NIR'])

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
    p.darkCorrect()
    p.flatCorrect()
    p.separateSky()
    p.associateSky()
    p.skyCorrect(mask_objects=False)
    p.alignAndStack()
    return

default = reduce
