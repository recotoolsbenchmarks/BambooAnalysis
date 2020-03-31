"""
The :py:mod:`bamboo.plots` module provides high-level classes to represent
and manipulate selections and plots.
"""
__all__ = ("Plot", "EquidistantBinning", "VariableBinning", "Selection", "DerivedPlot", "SummedPlot")

import logging
logger = logging.getLogger(__name__)

from itertools import chain
from . import treeoperations as top

class FactoryBackend(object):
    """ Interface for factory backend (to separate Plots classes from ROOT::RDataFrame part) """
    def __init__(self):
        pass
    def addSelection(self, selection):
        pass
    def addPlot(self, plot):
        pass
    def addDerived(self, product):
        pass
    def getResults(self, plot):
        pass
    def writeSkim(self, sele, outputFile, treeName, definedBranches=None, origBranchesToKeep=None, maxSelected=-1):
        pass
    @staticmethod
    def create(tree):
        """ Factory method, should return a pair of the backend and root selection """
        return (None, None)

class Product:
    """ Interface for output products (plots, counters etc.) """
    __slots__ = ("name",)
    def __init__(self, name):
        self.name = name
    def produceResults(self, bareResults, fbe):
        """
        Main interface method, called by the backend

        :param bareResults: iterable of histograms for this plot produced by the backend
        :param fbe: reference to the backend

        :returns: an iterable with ROOT objects to save to the output file
        """
        pass

class EquidistantBinning:
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

class VariableBinning:
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

class Plot(Product):
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
    def __init__(self, name, variables, selection, binnings, title="", axisTitles=tuple(), axisBinLabels=tuple(), plotopts=None, autoSyst=True):
        """ Generic constructor. Please use the static :py:meth:`~bamboo.plots.Plot.make1D`,
        :py:meth:`~bamboo.plots.Plot.make2D` and :py:meth:`~bamboo.plots.Plot.make3D` methods,
        which provide a more convenient interface to construct histograms
        (filling in some defaults requires knowing the dimensionality).
        """
        if len(variables) != len(binnings):
            raise ValueError("Unequal number of variables ({0:d}) and binnings ({1:d})".format(len(variables), len(binnings)))
        super(Plot, self).__init__(name)
        self.variables = variables
        self.selection = selection
        self.binnings = binnings
        self.title = title
        self.axisTitles = axisTitles
        self.axisBinLabels = axisBinLabels
        self.plotopts = plotopts if plotopts else dict()
        ## register with backend
        selection._fbe.addPlot(self, autoSyst=autoSyst)

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

    def produceResults(self, bareResults, fbe):
        """
        Trivial implementation of :py:meth:`~bamboo.plots.Product.produceResults`, return ``bareResults``

        Subclasses can e.g. calculate additional systematic variation histograms from the existing ones

        :param bareResults: list of nominal and systematic variation histograms for this :py:class:`~bamboo.plots.Plot`
        :param fbe: reference to the backend

        :returns: ``bareResults``
        """
        return bareResults

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
        :param xTitle: x-axis title (optional, taken from plot title by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)
        :param autoSyst: automatically add systematic variations (True by default - set to False to turn off)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 1-dimensional histogram

        :Example:

        >>> hasTwoEl = noSel.refine(cut=(op.rng_len(t.Electron) >= 2))
        >>> mElElPlot = Plot.make1D("mElEl", op.invariant_mass(t.Electron[0].p4, t.Electron[1].p4), hasTwoEl, EquidistantBinning(80, 50., 130.), title="Invariant mass of the leading-PT electrons")
        """
        title = kwargs.pop("title", "")
        kwargs["axisTitles"] = (kwargs.pop("xTitle", title),)
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None),)
        return Plot(name, (adaptArg(variable),), selection, (binning,), **kwargs)
    @staticmethod
    def make2D(name, variables, selection, binnings, **kwargs):
        """ Construct a 2-dimensional histogram plot

        :param name: unique plot name
        :param variables: x- and y-axis variable expression (iterable, e.g. tuple or list)
        :param selection: the :py:class:`~bamboo.plots.Selection` with cuts and weights to apply
        :param binnings: x- and y-axis binnings (iterable, e.g. tuple or list)
        :param title: plot title
        :param xTitle: x-axis title (optional, empty by default)
        :param yTitle: y-axis title (optional, empty by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param yBinLabels: y-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)
        :param autoSyst: automatically add systematic variations (True by default - set to False to turn off)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 2-dimensional histogram
        """
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)
    @staticmethod
    def make3D(name, variables, selection, binnings, **kwargs):
        """ Construct a 3-dimensional histogram

        :param name: unique plot name
        :param variables: x-, y- and z-axis variable expression (iterable, e.g. tuple or list)
        :param selection: the :py:class:`~bamboo.plots.Selection` with cuts and weights to apply
        :param binnings: x-, y-, and z-axis binnings (iterable, e.g. tuple or list)
        :param title: plot title
        :param xTitle: x-axis title (optional, empty by default)
        :param yTitle: y-axis title (optional, empty by default)
        :param zTitle: z-axis title (optional, empty by default)
        :param xBinLabels: x-axis bin labels (optional)
        :param yBinLabels: y-axis bin labels (optional)
        :param zBinLabels: z-axis bin labels (optional)
        :param plotopts: dictionary of options to pass directly to plotIt (optional)
        :param autoSyst: automatically add systematic variations (True by default - set to False to turn off)

        :returns: the new :py:class:`~bamboo.plots.Plot` instance with a 3-dimensional histogram
        """
        kwargs["axisTitles"] = (kwargs.pop("xTitle", ""), kwargs.pop("yTitle", ""), kwargs.pop("zTitle", ""))
        kwargs["axisBinLabels"] = (kwargs.pop("xBinLabels", None), kwargs.pop("yBinLabels", None), kwargs.pop("zBinLabels", None))
        return Plot(name, tuple(adaptArg(v) for v in variables), selection, binnings, **kwargs)

class Selection:
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
    def __init__(self, parent, name, cuts=None, weights=None, autoSyst=True):
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
        self._cSysts = top.collectSystVars(self._cuts)
        self._wSysts = top.collectSystVars(self._weights)

        ## register with backend
        if isinstance(parent, Selection):
            self.autoSyst = parent.autoSyst and autoSyst
            self.parent = parent
            self._fbe = parent._fbe
        else:
            self.autoSyst = autoSyst
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
    @property
    def weightSystematics(self):
        if self.parent:
            return top.mergeSystVars(self.parent.weightSystematics, self._wSysts)
        else:
            return dict(self._wSysts)
    @property
    def cutSystematics(self):
        if self.parent:
            return top.mergeSystVars(self.parent.cutSystematics, self._cSysts)
        else:
            return dict(self._cSysts)
    @property
    def systematics(self):
        return top.mergeSystVars(self.weightSystematics, self.cutSystematics)
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
        return ( ( len(self.cuts) == len(other.cuts) ) and all( sc == oc for sc,oc in zip(self.cuts, other.cuts) )
             and ( len(self.weights) == len(other.weights) ) and all( sw == ow for sw,ow in zip(self.weights, other.weights) ) )

    def refine(self, name, cut=None, weight=None, autoSyst=True):
        """ Create a new selection by adding cuts and/or weight factors

        :param name: unique name of the new selection
        :param cut: expression (or list of expressions) with additional selection criteria (combined using logical AND)
        :param weight: expression (or list of expressions) with additional weight factors
        :param autoSyst: automatically add systematic variations (True by default - set to False to turn off; note that this would also turn off automatic systematic variations for any selections and plots that derive from the one created by this method)

        :returns: the new :py:class:`~bamboo.plots.Selection`
        """
        return Selection(self, name,
                cuts   =( ( adaptArg(ct, "Bool_t") for ct in (cut    if hasattr(cut   , "__len__") else [cut   ]) ) if cut    else None ),
                weights=( ( adaptArg(wt, "Bool_t") for wt in (weight if hasattr(weight, "__len__") else [weight]) ) if weight else None ),
                autoSyst=autoSyst
                )

    @staticmethod
    def _makeExprAnd(listOfReqs):
        # op.AND for expressions (helper for histfactory etc.)
        if len(listOfReqs) > 1:
            return adaptArg(op.AND(*listOfReqs))
        elif len(listOfReqs) == 1:
            return listOfReqs[0]
        else:
            from .treeproxies import makeConst, boolType
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

class DerivedPlot(Product):
    """
    Base class for a plot with results based on other plots' results

    The :py:attr:`~bamboo.plots.DerivedPlot.dependencies` attribute that lists
    the :py:class:`~bamboo.plots.Plot`-like objects this one depends on (which
    may be used e.g. to order operations).
    The other necessary properties (binnings, titles, labels, etc.) are taken
    from the keyword arguments to the constructor, or the first dependency.
    The :py:meth:`~bamboo.plots.DerivedPlot.produceResults` method,
    which is called by the backend to retrieve the derived results,
    should be overridden with the desired calculation.

    Typical use cases are summed histograms, background subtraction, etc.
    (the results are combined for different subjobs with hadd, so derived 
    quantities that require the full statistics should be calculated from
    the postprocessing step; alternative or additional systematic variations
    calculated from the existing ones can be added by subclassing
    :py:class:`~bamboo.plots.Plot`).
    """
    def __init__(self, name, dependencies, **kwargs):
        super(DerivedPlot, self).__init__(name)
        self.dependencies = dependencies
        self.binnings = kwargs.get("binnings", dependencies[0].binnings)
        self.axisTitles = kwargs.get("axisTitles",
                tuple([ kwargs.get("{0}Title".format(ax), dependencies[0].axisTitles[i])
                    for i,ax in enumerate("xyzuvw"[:len(self.variables)]) ]))
        self.axisBinLabels = kwargs.get("axisBinLabels",
                tuple([ kwargs.get("{0}BinLabels".format(ax), dependencies[0].axisBinLabels[i])
                    for i,ax in enumerate("xyzuvw"[:len(self.variables)]) ]))
        self.plotopts = kwargs.get("plotopts", dict())
        # register with backend
        dependencies[0].selection._fbe.addDerived(self)
    @property
    def variables(self):
        return [ None for x in self.binnings ]
    def produceResults(self, bareResults, fbe):
        """
        Main interface method, called by the backend

        :param bareResults: iterable of histograms for this plot produced by the backend (none)
        :param fbe: reference to the backend, can be used to retrieve the histograms for the dependencies, e.g. with :py:meth:`~bamboo.plots.DerivedPlot.collectDependencyResults`

        :returns: an iterable with ROOT objects to save to the output file
        """
        return []
    def collectDependencyResults(self, fbe):
        """ helper method: collect all results of the dependencies

        :returns: ``[ (nominalResult, {"variation" : variationResult}) ]``
        """
        res_dep = []
        for dep in self.dependencies:
            resNom = None
            resPerVar = {}
            for res in fbe.getResults(dep):
                if "__" not in res.GetName():
                    assert resNom is None
                    resNom = res
                else:
                    resVar = res.GetName().split("__")[1]
                    resPerVar[resVar] = res
            res_dep.append((resNom, resPerVar))
        return res_dep

class SummedPlot(DerivedPlot):
    """ A :py:class:`~bamboo.plots.DerivedPlot` implementation that sums histograms """
    def __init__(self, name, termPlots, **kwargs):
        super(SummedPlot, self).__init__(name, termPlots, **kwargs)
    def produceResults(self, bareResults, fbe):
        res_dep = self.collectDependencyResults(fbe)
        # list all variations (some may not be there for all)
        allVars = set()
        for _,resVar in res_dep:
            allVars.update(resVar.keys())
        # sum nominal
        hNom = res_dep[0][0].Clone(self.name)
        for ihn,_ in res_dep[1:]:
            hNom.Add(ihn.GetPtr())
        results = [ hNom ]
        # sum variations (using nominal if not present for some)
        for vn in allVars:
            hVar = res_dep[0][1].get(vn, res_dep[0][0]).Clone("__".join((self.name, vn)))
            for ihn,ihv in res_dep[1:]:
                hVar.Add(ihv.get(vn, ihn).GetPtr())
            results.append(hVar)
        return results

class CutFlowReport(Product):
    class Entry:
        def __init__(self, name, nominal=None, systVars=None, parent=None, children=None):
            self.name = name
            self.nominal = nominal
            self.systVars = systVars
            self.parent = parent
            self.children = list(children) if children is not None else []
        def setParent(self, parent):
            self.parent = parent
            if self not in parent.children:
                parent.children.append(self)
    def __init__(self, name, selections, recursive=True, autoSyst=False, cfres=None):
        super(CutFlowReport, self).__init__(name)
        self.recursive = recursive
        self.selections = list(selections) if hasattr(selections, "__iter__") else [selections]
        self.cfres = cfres
        if selections and cfres is None:
            self.selections[0]._fbe.addCutFlowReport(self, autoSyst=autoSyst)
    def produceResults(self, bareResults, fbe):
        self.cfres = fbe.results[self.name]
        entries = list()
        for stat in self.cfres:
            iEn = stat
            while iEn is not None and iEn not in entries:
                entries.append(iEn)
                iEn = iEn.parent
        return [ res.nominal for res in entries ] + list(chain.from_iterable(res.systVars.values() for res in entries))
    def rootEntries(self):
        ## helper: traverse reports tree up
        def travUp(entry):
            yield entry
            yield from travUp(entry.parent)
        return set(next(en for en in travUp(res) if en.parent is None) for res in self.cfres)

    def readFromResults(self, resultsFile):
        cfres = []
        entries = dict() # by selection name
        for sel in self.selections:
            if sel.name not in entries:
                entries[sel.name] = CutFlowReport.Entry(sel.name)
                if self.recursive:
                    isel = sel.parent
                    entry_d = entries[sel.name]
                    while isel is not None:
                       if isel.name in entries:
                           entry_p = entries[isel.name]
                           entry_d.setParent(entry_p)
                           break
                       entry_p = CutFlowReport.Entry(isel.name)
                       entries[isel.name] = entry_p
                       entry_d.setParent(entry_p)
                       entry_d = entry_p
                       isel = isel.parent
            cfres.append(entries[sel.name])
        ## retrieve nominal
        toRem = []
        for selName, entry in entries.items():
            kyNm = f"{self.name}_{selName}"
            entry.nominal = resultsFile.Get(kyNm)
            if not entry.nominal:
                logger.warning(f"Could not find object {kyNm}")
                toRem.append(selName)
        ## remove if not found
        for ky in toRem:
            if entries[ky].parent:
                for ch in entries[ky].children:
                    ch.setParent(entries[ky].parent)
                entries[ky].parent.children.remove(entries[ky])
            del entries[ky]
        ## and systematic variations
        prefix = f"{self.name}_"
        for ky in resultsFile.GetListOfKeys():
            if ky.GetName().startswith(prefix):
                selName = ky.GetName().split("__")[0][len(prefix):]
                if selName in entries:
                    entry = entries[selName]
                    cnt = ky.GetName().count("__")
                    if cnt == 1:
                        varNm = ky.GetName().split("__")[1]
                        if varNm in entry.systVars:
                            logger.warning(f"{self.name}: counter for variation {varNm} already present for selection {selName}")
                        entry.systVars[varNm] = ky.ReadObject()
                    elif cnt > 1:
                        logger.warning("Key {ky.GetName()!r} contains '__' more than once, this will break assumptions")
        return CutFlowReport(self.name, self.selections, recursive=self.recursive, cfres=cfres)
