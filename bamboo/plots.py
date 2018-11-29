"""
Helper classes describing histograms, binnings etc.
"""
__all__ = ("Plot", "EquidistantBinning", "VariableBinning", "Selection")

import numpy as np
from itertools import chain

class FactoryBackend(object):
    """ Interface for factory backend (to separate Plots classes from ROOT::RDataFrame part) """
    def __init__(self):
        pass
    def addSelection(self, selection):
        pass
    def addPlot(self, plot):
        pass
    @staticmethod
    def create(tree):
        """ Factory method, should return a pair of the backend and root selection """
        return (None, None)

class EquidistantBinning(object):
    """ Equidistant binning
    """
    __slots__ = ("__weakref__", "N", "mn", "mx")
    def __init__(self, N, mn, mx):
        self.N = N
        self.mn = mn
        self.mx = mx
    @property
    def minimum(self):
        return self.mn
    @property
    def maximum(self):
        return self.mx

class VariableBinning(object):
    """ Variable-sized bins
    """
    __slots__ = ("__weakref__", "binEdges")
    def __init__(self, binEdges):
        self.binEdges = np.array(binEdges)
    @property
    def N(self):
        return len(self.binEdges)-1
    @property
    def minimum(self):
        return self.binEdges[0]
    @property
    def maximum(self):
        return self.binEdges[-1]

from .treeproxies import adaptArg
from . import treefunctions as op

class Plot(object):
    """ Plot representation (specifications for making a 1,2,3,...-dimensional histogram) """
    __slots__ = ("__weakref__", "name", "variables", "selection", "binnings", "title", "axisTitles", "axisBinLabels", "plotopts", "df")
    def __init__(self, name, variables, selection, binnings, title="", axisTitles=tuple(), axisBinLabels=tuple(), plotopts=None):
        if len(variables) != len(binnings):
            raise ValueError("Unequal number of variables ({0:d}) and binnings ({1:d})".format(len(variables), len(binnings)))
        self.name = name
        self.variables = variables
        self.selection = selection
        self.binnings = binnings
        self.title = title
        self.axisTitles = axisTitles
        self.axisBinLabels = axisBinLabels
        self.plotopts = plotopts if plotopts else dict()
        ## register with backend
        selection._fbe.addPlot(self)

    def clone(self, name=None, variables=None, selection=None, binnings=None, title=None, axisTitles=None, axisBinLabels=None, plotopts=None):
        """ Low-level helper method: copy with optional re-setting of attributes
        """
        return Plot( (name if name is not None else self.name)
                   , (variables if variables is not None else self.variables)
                   , (selection if selection is not None else self.selection)
                   , (binnings if binnings is not None else self.binnings)
                   , title=(title if title is not None else self.title)
                   , axisTitles=(axisTitles if axisTitles is not None else self.axisTitles)
                   , axisBinLabels=(axisBinLabels if axisBinLabels is not None else self.axisBinLabels)
                   , plotopts=(plotopts if plotopts is not None else self.plotopts)
                   )

    @property
    def cut(self):
        return self.selection.cut
    @property
    def weight(self):
        return self.selection.weight

    @property
    def longTitle(self):
        return ";".join(chain([self.title], self.axisTitles))

    def __repr__(self):
        return "Plot({0!r}, {1!r}, {2!r}, {3!r}, title={4!r}, axisTitles={5!r})".format(self.name, self.variables, self.selection, self.binnings, self.title, self.axisTitles)

    @staticmethod
    def make1D(name, variable, selection, binning, **kwargs):
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""),)
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None),)
        return Plot(name, (adaptArg(variable),), selection, (binning,), **kwargs)
    @staticmethod
    def make2D(name, variables, selection, binnings, **kwargs):
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)
    @staticmethod
    def make3D(name, variables, selection, binnings, **kwargs):
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""), kwargs.pop("zTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None), kwargs.pop("zBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)

class Selection(object):
    """ Group of selection requirements and weight factors
    """
    def __init__(self, parent, name, cuts=None, weights=None):
        self.name      = name
        self.parent   = None
        self._cuts     = [ adaptArg(cut, "Bool_t") for cut in cuts ] if cuts else []
        self._weights  = [ adaptArg(wgt, typeHint="Float_t") for wgt in weights ] if weights else []
        ## register with backend
        if isinstance(parent, Selection):
            self.parent = parent
            self._fbe = parent._fbe
        else:
            assert isinstance(parent, FactoryBackend)
            self._fbe = parent
        self._fbe.addSelection(self)

    @property
    def cuts(self):
        if self.parent:
            return self.parent.cuts + self._cuts
        else:
            return self._cuts
    @property
    def weights(self):
        if self.parent:
            return self.parent.weights + self._weights
        else:
            return self._weights
    ## for debugging/monitoring: full cut and weight expression ## TODO review
    @property
    def cut(self):
        return Selection._makeExprAnd(self.cuts)
    @property
    def weight(self):
        return Selection._makeExprProduct(self.weights)
    def __repr__(self): ## TODO maybe change to refer to parent
        return "Selection({0!r}, {1!r}, {2!r})".format(self.cut, self.weight)
    def __eq__(self, other):
        ## FIXME do we even still need this?
        return ( ( len(self.cuts) == len(other.cuts) ) and all( sc == oc for sc,oc in izip(self.cuts, other.cuts) )
             and ( len(self.weights) == len(other.weights) ) and all( sw == ow for sw,ow in izip(self.weights, other.weights) ) )

    def refine(self, name, cut=None, weight=None):
        """ Create a new selection by adding a cut and/or weight """
        return Selection(self, name,
                cuts   =( ( adaptArg(ct, "Bool_t") for ct in (cut    if hasattr(cut   , "__len__") else [cut   ]) ) if cut    else None ),
                weights=( ( adaptArg(wt, "Bool_t") for wt in (weight if hasattr(weight, "__len__") else [weight]) ) if weight else None )
                )

    @staticmethod
    def _makeExprAnd(listOfReqs):
        """ op.AND for expressions (helper for histfactory etc.)
        """
        if len(listOfReqs) > 1:
            return adaptArg(op.AND(*listOfReqs))
        elif len(listOfReqs) == 1:
            return listOfReqs[0]
        else:
            from .treeproxies import makeConst
            return adaptArg(makeConst("true", boolType), typeHint=boolType)
    @staticmethod
    def _makeExprProduct(listOfFactors):
        """ op.product for expressions (helper for histfactory etc.)
        """
        if len(listOfFactors) > 1:
            return adaptArg(op.product(*listOfFactors))
        elif len(listOfFactors) == 1:
            return listOfFactors[0]
        else:
            from .treeproxies import makeConst
            return adaptArg(makeConst(1., "float"), typeHint="float")
