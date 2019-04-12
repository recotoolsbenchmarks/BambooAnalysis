"""
ROOT::RDataFrame backend classes
"""
import logging
logger = logging.getLogger(__name__)

from itertools import chain, count

from .plots import FactoryBackend, Selection
from . import treefunctions as op
from . import treeoperations as top

class SelWithDefines(top.CppStrRedir):
    def __init__(self, parent, df, weights=None, wName=None):
        self.df = df
        self.explDefine = list()
        if isinstance(parent, SelWithDefines):
            self.parent = parent
            self.backend = parent.backend
            self._definedColumns = dict(parent._definedColumns)
        else:
            self.parent = None
            self.backend = parent
            self._definedColumns = dict()

        if weights:
            self._initWeights(wName, weights)
            self.wName = wName
        elif self.parent and self.parent.wName:
            self.wName = self.parent.wName
        else:
            assert not wName
            self.wName = None
        top.CppStrRedir.__init__(self)
    
    def _initWeights(self, wName, weights):
        weightExpr = Selection._makeExprProduct(
            ([top.adaptArg(op.extVar("float", self.parent.wName), typeHint="float")]+weights) if self.parent.wName
            else weights
            )
        self._define(wName, weightExpr)

    def _getColName(self, op):
        if op in self._definedColumns:
            return self._definedColumns[op]

    def define(self, expr):
        """ explicitly define column for expression (returns the column name) """
        if not self._getColName(expr):
            self.explDefine.append(expr)
        return self(expr)

    def _define(self, name, expr):
        logger.debug("Defining {0} as {1}".format(name, expr.get_cppStr(defCache=self)))
        self.df = self.df.Define(name, expr.get_cppStr(defCache=self))
        self._definedColumns[expr] = name

    def symbol(self, decl, **kwargs):
        return self.backend.symbol(decl, **kwargs)

    def _inExplDefines(self, arg):
        return arg in self.explDefine or ( self.parent and self.parent._inExplDefines(arg) )

    def shouldDefine(self, arg):
        return self.backend.shouldDefine(arg, defCache=self) or self._inExplDefines(arg)

    def __call__(self, arg):
        """ Get the C++ string corresponding to an op """
        if not self.shouldDefine(arg): ## the easy case
            try:
                return arg.get_cppStr(defCache=self)
            except Exception as ex:
                logger.error("Could not get cpp string for {0!r}: {1!r}".format(arg, ex))
                return "NONE"
        else:
            nm = self._getColName(arg)
            if not nm: ## define it then
                nm = self.backend.getUColName()
                self._define(nm, arg)
            return nm

## NOTE these are global (for the current process&interpreter)
## otherwise they would be overwritten in sequential mode (this even allows reuse)
_gSymbols = dict()
_giFun = 0

class DataframeBackend(FactoryBackend):
    def __init__(self, tree, outFileName=None):
        from cppyy import gbl
        self.rootDF = gbl.RDataFrame(tree)
        self.outFile = gbl.TFile.Open(outFileName, "CREATE") if outFileName else None
        self.selDFs = dict()      ## selection name -> SelWithDefines
        self.plotResults = dict() ## plot name -> result pointer
        super(DataframeBackend, self).__init__()
        self._iCol = 0
    def _getUSymbName(self):
        global _giFun
        _giFun += 1
        return "myFun{0:d}".format(_giFun)
    def getUColName(self):
        self._iCol += 1
        return "myCol{0:d}".format(self._iCol)

    def shouldDefine(self, op, defCache=None):
        return ( any(isinstance(op, expType) for expType in (top.Select, top.Sort, top.Map, top.Next, top.Reduce, top.Combine))
                and not any(op.deps(defCache=defCache, select=lambda dp : isinstance(dp, top.LocalVariablePlaceholder))) )

    def symbol(self, decl, resultType=None, args=None, nameHint=None):
        global _gSymbols
        if decl in _gSymbols:
            return _gSymbols[decl]
        else:
            if nameHint and nameHint not in _gSymbols.values():
                name = nameHint
            else:
                name = self._getUSymbName()
            _gSymbols[decl] = name
            if resultType and args: ## then it needs to be wrapped in a function
                fullDecl = "{result} {name}({args})\n{{\n  return {0};\n}};\n".format(
                            decl, result=resultType, name=name, args=args)
            else:
                fullDecl = decl.replace("<<name>>", name)

            logger.debug("Defining new symbol with interpreter: {0}".format(fullDecl))
            from cppyy import gbl
            gbl.gInterpreter.Declare(fullDecl)
            return name

    @staticmethod
    def create(decoTree, outFileName=None):
        inst = DataframeBackend(decoTree._tree, outFileName=None)
        rootSel = Selection(inst, "none")
        return inst, rootSel

    def addSelection(self, sele):
        """ Define ROOT::RDataFrame objects needed for this selection """
        if sele.name in self.selDFs:
            raise ValueError("A Selection with the name '{0}' already exists".format(sele.name))
        parentDF = self.selDFs[sele.parent.name] if sele.parent else None
        if sele._cuts:
            assert parentDF ## FIXME if not something went wrong - there *needs* to be a root no-op sel
            expr = Selection._makeExprAnd(sele._cuts)
            filterStr = expr.get_cppStr(defCache=parentDF)
            logger.debug("Filtering with {0}".format(filterStr))
            selDF = parentDF.df.Filter(filterStr)
        else:
            if parentDF:
                selDF = parentDF.df
            else:
                selDF = self.rootDF

        self.selDFs[sele.name] = SelWithDefines((parentDF if sele.parent else self), selDF, weights=sele._weights, wName=("w_{0}".format(sele.name) if sele._weights else None))

    def addPlot(self, plot):
        """ Define ROOT::RDataFrame objects needed for this plot (and keep track of the result pointer) """
        if plot.name in self.plotResults:
            raise ValueError("A Plot with the name '{0}' already exists".format(plot.name))
        ## TODO DataFrame might throw (or segfault), but we should catch all possible errors before
        ## NOTE pre/postfixes should already be inside the plot name
        hModel = DataframeBackend.makePlotModel(plot)
        selND = self.selDFs[plot.selection.name]
        varExprs = dict(("v{0:d}_{1}".format(i, plot.name), selND(var)) for i,var in zip(count(), plot.variables))
        plotDF = selND.df ## after getting the expressions, to pick up columns that were defined on-demand
        for vName, vExpr in varExprs.items():
            #logger.debug("Defining {0} as {1} (defined for {2}: {3})".format(vName, vExpr, plotDF, ", ".join(plotDF.GetDefinedColumnNames()))) ## needs 6.16
            logger.debug("Defining {0} as {1}".format(vName, vExpr))
            plotDF = plotDF.Define(vName, vExpr)
        plotFun = getattr(plotDF, "Histo{0:d}D".format(len(plot.variables)))
        if selND.wName:
            logger.debug("Adding plot {0} with variables {1} and weight {2}".format(plot.name, ", ".join(varExprs.keys()), selND.wName))
            plotDF = plotFun(hModel, *chain(varExprs.keys(), [selND.wName]))
        else:
            logger.debug("Adding plot {0} with variables {1}".format(plot.name, ", ".join(varExprs.keys())))
            plotDF = plotFun(hModel, *varExprs.keys())
        self.plotResults[plot.name] = plotDF

    @staticmethod
    def makePlotModel(plot):
        from cppyy import gbl
        modCls = getattr(gbl.RDF, "TH{0:d}DModel".format(len(plot.binnings)))
        return modCls(plot.name, plot.title, *chain.from_iterable(
            DataframeBackend.makeBinArgs(binning) for binning in plot.binnings))
    @staticmethod
    def makeBinArgs(binning):
        from .plots import EquidistantBinning, VariableBinning
        if isinstance(binning, EquidistantBinning):
            return (binning.N, binning.mn, binning.mx)
        elif isinstance(binning, VariableBinning):
            return (binning.N, binning.binEdges)
        else:
            raise ValueError("Binning of unsupported type: {0!r}".format(binning))

    def getPlotResult(self, plot):
        return self.plotResults[plot.name]
