"""
The :py:mod:`bamboo.plots` module provides high-level classes to represent
and manipulate selections and plots.
"""
__all__ = ("Plot", "EquidistantBinning", "VariableBinning", "Selection")

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
    """ Equidistant binning """
    __slots__ = ("__weakref__", "N", "mn", "mx")
    def __init__(self, N, mn, mx):
        """
        :param N: number of bins
        :param mn: minimum axis value
        :param mx: maximum axis value
        """
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
    """ Variable-sized binning """
    __slots__ = ("__weakref__", "binEdges")
    def __init__(self, binEdges):
        """
        :param binEdges: iterable with the edges. There will be ``len(binEges)-1`` bins
        """
        self.binEdges = list(binEdges)
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
    """ A :py:class:`~bamboo.plots.Plot` object contains all information needed
    to produce a histogram: the variable(s) to plot, binnings and options
    (axis titles, optionally some style information), and a reference to
    a :py:class:`~bamboo.plots.Selection` (which holds all cuts and weights to apply for the plot).

    .. note::

        All :py:class:`~bamboo.plots.Plot` (and :py:class:`~bamboo.plots.Selection`) instances
        need to have a unique name. This name is used to construct output filenames, and internally
        to define DataFrame columns with readable names.
        The constructor will raise an exception if an existing name is used.
    """
    __slots__ = ("__weakref__", "name", "variables", "selection", "binnings", "title", "axisTitles", "axisBinLabels", "plotopts", "df")
    def __init__(self, name, variables, selection, binnings, title="", axisTitles=tuple(), axisBinLabels=tuple(), plotopts=None):
        """ Generic constructor. Please use the static :py:meth:`~bamboo.plots.Plot.make1D`,
        :py:meth:`~bamboo.plots.Plot.make2D` and :py:meth:`~bamboo.plots.Plot.make3D` methods,
        which provide a more convenient interface to construct histograms
        (filling in some defaults requires knowing the dimensionality).
        """
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
        """ Helper method: create a copy with optional re-setting of attributes """
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
        """ Construct a 1-dimensional histogram plot

        :param name: unique plot name
        :param variable: x-axis variable expression
        :param selection: the :py:class:`~bamboo.plots.Selection` with cuts and weights to apply
        :param binning: x-axis binning
        :param title: plot title
        :param xTitle: x-axis title (optional, empty by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 1-dimensional histogram

        :Example:

        >>> hasTwoEl = noSel.refine(cut=(op.rng_len(t.Electron) >= 2))
        >>> mElElPlot = Plot.make1D("mElEl", op.invariant_mass(t.Electron[0].p4, t.Electron[1].p4), hasTwoEl, EquidistantBinning(80, 50., 130.), title="Invariant mass of the leading-PT electrons")
        """
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""),)
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None),)
        return Plot(name, (adaptArg(variable),), selection, (binning,), **kwargs)
    @staticmethod
    def make2D(name, variables, selection, binnings, **kwargs):
        """ Construct a 2-dimensional histogram plot

        :param name: unique plot name
        :param variable: x-axis variable expression
        :param selection: the :py:class:`~bamboo.plots.Selection` with cuts and weights to apply
        :param binnings: x- and y-axis binnings
        :param title: plot title
        :param xTitle: x-axis title (optional, empty by default)
        :param yTitle: y-axis title (optional, empty by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param yBinLabels: y-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 2-dimensional histogram
        """
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)
    @staticmethod
    def make3D(name, variables, selection, binnings, **kwargs):
        """ Construct a 3-dimensional histogram

        :param name: unique plot name
        :param variable: x-axis variable expression
        :param selection: the :py:class:`~bamboo.plots.Selection` with cuts and weights to apply
        :param binnings: x-, y-, and z-axis binnings
        :param title: plot title
        :param xTitle: x-axis title (optional, empty by default)
        :param yTitle: y-axis title (optional, empty by default)
        :param zTitle: z-axis title (optional, empty by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param yBinLabels: y-axis bin labels (optional)
        :param zBinLabels: z-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 3-dimensional histogram
        """
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""), kwargs.pop("zTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None), kwargs.pop("zBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)

class Selection(object):
    """ A :py:class:`~bamboo.plots.Selection` object groups a set of selection criteria
    (cuts) and weight factors that belong to a specific stage of the selection and analysis.
    Selections should be constructed by calling the :py:meth:`~bamboo.plots.Selection.refine`
    method on a 'root' selection (which may include overall selections and weights, e.g.
    a lumi mask for data and pileup reweighting for MC).

    .. note::

        All :py:class:`~bamboo.plots.Selection` (and :py:class:`~bamboo.plots.Plot`) instances
        need to have a unique name. This name is used internally to define DataFrame columns
        with readable names.
        The constructor will raise an exception if an existing name is used.
    """
    def __init__(self, parent, name, cuts=None, weights=None):
        """ Constructor. Prefer using :py:meth:`~bamboo.plots.Selection.refine` instead (except for the 'root' selection)

        :param parent: backend or parent selection
        :param name: (unique) name of the selection
        :param cuts: iterable of selection criterion expressions (optional)
        :param weights: iterable of weight factors (optional)
        """
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
        return "Selection({0!r}, {1!r}, {2!r})".format(self.name, self.cut, self.weight)
    def __eq__(self, other):
        ## FIXME do we even still need this?
        return ( ( len(self.cuts) == len(other.cuts) ) and all( sc == oc for sc,oc in izip(self.cuts, other.cuts) )
             and ( len(self.weights) == len(other.weights) ) and all( sw == ow for sw,ow in izip(self.weights, other.weights) ) )

    def refine(self, name, cut=None, weight=None):
        """ Create a new selection by adding a cuts and/or weight factors

        :param name: unique name of the new selection
        :param cut: expression (or list of expressions) with additional selection criteria
        :param weight: expression (or list of expressions) with additional weight factors

        :returns: the new :py:class:`~bamboo.plots.Selection`
        """
        return Selection(self, name,
                cuts   =( ( adaptArg(ct, "Bool_t") for ct in (cut    if hasattr(cut   , "__len__") else [cut   ]) ) if cut    else None ),
                weights=( ( adaptArg(wt, "Bool_t") for wt in (weight if hasattr(weight, "__len__") else [weight]) ) if weight else None )
                )

    @staticmethod
    def _makeExprAnd(listOfReqs):
        # op.AND for expressions (helper for histfactory etc.)
        if len(listOfReqs) > 1:
            return adaptArg(op.AND(*listOfReqs))
        elif len(listOfReqs) == 1:
            return listOfReqs[0]
        else:
            from .treeproxies import makeConst
            return adaptArg(makeConst("true", boolType), typeHint=boolType)
    @staticmethod
    def _makeExprProduct(listOfFactors):
        # op.product for expressions (helper for histfactory etc.)
        if len(listOfFactors) > 1:
            return adaptArg(op.product(*listOfFactors))
        elif len(listOfFactors) == 1:
            return listOfFactors[0]
        else:
            from .treeproxies import makeConst
            return adaptArg(makeConst(1., "float"), typeHint="float")
