"""
Recipes available to data with tags ['OCTOCAM', 'SPECT', 'NIR'].
Default is "reduce".
"""
recipe_tags = set(['OCTOCAM', 'SPECT', 'NIR'])

def reduce(p):
    """
    This recipe process near-IR spectrum up to and including alignment and
    stacking.  A single stacked extracted and calibrated spectrum is produced.


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
    # whatever else is needed.
    # ...
    return

default = reduce
