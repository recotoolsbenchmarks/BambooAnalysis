"""
Analysis helper functions that don't fit elsewhere
"""
from . import treefunctions as op

def addLumiMask(sel, jsonName, runRange=None, runAndLS=None, name="goodlumis"):
    """ Refine selection with a luminosity block filter

    Typically applied directly to the root selection (for data).
    runAndLS should be a tuple of expressions with the run number and luminosity block ID.
    The run range is used to limit the part of the JSON file to consider,
    see the LumiMask helper class for details.
    """
    if runAndLS is None:
        raise RuntimeError("Cannot construct a filter for the good lumi blocks without accessors (backend.create(..., runAndLS=XXX)), tree->(run, LS)")
    lumiSel = op.define("LumiMask", 'const auto <<name>> = LumiMask::fromJSON("{0}"{1});'.format(
                jsonName, (", {0:d}, {1:d}".format(*runRange) if runRange is not None else "")))
    return sel.refine(name, cut=lumiSel.accept(*runAndLS))
