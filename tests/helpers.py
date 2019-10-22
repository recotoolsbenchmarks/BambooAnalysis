import math

def isclose_float(a, b):
    from bamboo.root import gbl
    return math.isclose(a, b, rel_tol=getattr(gbl, "std::numeric_limits<float>").epsilon())
