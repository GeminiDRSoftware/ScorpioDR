"""
Recipes available to data with tags ['SCORPIO', 'BUNDLE'].
Default is "makeProcessedUnbundle".
"""
recipe_tags = set(['SCORPIO', 'BUNDLE', 'RAW', 'UNPREPARED'])

def makeProcessedUnbundle(p):
    """
    This recipe processes SCORPIO observation bundles.

    Parameters
    ----------
    p : PrimitivesCORE object
        A primitive set matching the recipe_tags.
    """

    p.separateChannels()
    return

_default = makeProcessedUnbundle
