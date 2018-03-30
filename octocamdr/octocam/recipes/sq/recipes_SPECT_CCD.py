"""
Recipes available to data with tags ['OCTOCAM', 'SPECT', 'CCD'].
Default is "reduce".
"""
recipe_tags = set(['OCTOCAM', 'SPECT', 'CCD'])

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
    p.flatCorrect()
    # whatever else is needed.
    # ...
    return

default = reduce
