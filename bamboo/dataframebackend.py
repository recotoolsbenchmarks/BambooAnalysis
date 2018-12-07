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
        self._definedColumns = dict()
        if isinstance(parent, SelWithDefines):
            self.parent = parent
            self.backend = self.parent.backend
        else:
            self.parent = None
            self.backend = parent
        if weights:
            self._initWeights(wName, weights)
        else:
            assert not wName
        self.wName = wName
        top.CppStrRedir.__init__(self)
    
    def _initWeights(self, wName, weights):
        weightExpr = Selection._makeExprProduct(
            ([top.adaptArg(op.extVar("float", self.parent.wName), typeHint="float")]+weights) if self.parent.wName
            else weights
            )
        logger.debug("Defining {0} as {1}".format(wName, self(weightExpr)))
        self.df = self.df.Define(wName, self(weightExpr))

    def _getColName(self, op):
        if op in self._definedColumns:
            return self._definedColumns[op]
        elif self.parent:
            res = self.parent._getColName(op)
            if res:
                return res

    def symbol(self, decl):
        return self.backend.symbol(decl)
    
    def __call__(self, arg):
        """ Get the C++ string corresponding to an op """
        if not self.backend.shouldDefine(arg): ## the easy case
            return arg.get_cppStr(defCache=self)
        else:
            nm = self._getColName(arg)
            if not nm: ## define it then
                nm = self.backend.getUColName()
                logger.debug("Defining {0} as {1}".format(nm, arg.get_cppStr(defCache=self)))
                self.df = self.df.Define(nm, arg.get_cppStr(defCache=self))
                self._definedColumns[arg] = nm
            return nm

class DataframeBackend(FactoryBackend):
    def __init__(self, tree, outFileName=None):
        import ROOT
        self.rootDF = ROOT.RDataFrame(tree)
        self.outFile = ROOT.TFile.Open(outFileName, "CREATE") if outFileName else None
        self.selDFs = dict()      ## selection name -> SelWithDefines
        self.plotResults = dict() ## plot name -> result pointer
        super(DataframeBackend, self).__init__()
        self._iCol = 0
        self._iFun = 0
        self._symbols = dict()
    def _getUSymbName(self):
        self._iFun += 1
        return "myFun{0:d}".format(self._iFun)
    def getUColName(self):
        self._iCol += 1
        return "myCol{0:d}".format(self._iCol)

    def shouldDefine(self, op):
        return any(isinstance(op, expType) for expType in (top.Select, top.Next, top.Reduce))

    def symbol(self, decl):
        if decl in self._symbols:
            return self._symbols[decl]
        else:
            name = self._getUSymbName()
            self._symbols[decl] = name
            fullDecl = decl.replace("<<name>>", name)
            logger.debug("Defining new symbol with interpreter: {0}".format(fullDecl))
            import ROOT
            ROOT.gInterpreter.Declare(fullDecl)
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
        if sele.parent:
            parentDF = self.selDFs[sele.parent.name]
            selDF = parentDF.df
        else:
            selDF = self.rootDF
        if sele._cuts:
            assert parentDF ## FIXME if not something went wrong - there *needs* to be a root no-op sel
            expr = Selection._makeExprAnd(sele._cuts)
            logger.debug("Filtering with {0}".format(expr.get_cppStr(defCache=parentDF)))
            selDF = selDF.Filter(expr.get_cppStr(defCache=parentDF))

        self.selDFs[sele.name] = SelWithDefines((parentDF if sele.parent else self), selDF, weights=sele._weights, wName=("w_{0}".format(sele.name) if sele._weights else None))

    def addPlot(self, plot):
        """ Define ROOT::RDataFrame objects needed for this plot (and keep track of the result pointer) """
        if plot.name in self.plotResults:
            raise ValueError("A Plot with the name '{0}' already exists".format(plot.name))
        ## TODO DataFrame might throw (or segfault), but we should catch all possible errors before
        ## NOTE pre/postfixes should already be inside the plot name
        hModel = DataframeBackend.makePlotModel(plot)
        selND = self.selDFs[plot.selection.name]
        plotDF = selND.df
        varNames = []
        for i,var in zip(count(), plot.variables):
            vName = "v{0:d}_{1}".format(i, plot.name)
            logger.debug("Defining {0} as {1}".format(vName, selND(var)))
            plotDF = plotDF.Define(vName, selND(var))
            varNames.append(vName)
        plotFun = getattr(plotDF, "Histo{0:d}D".format(len(plot.variables)))
        if selND.wName:
            plotDF = plotFun(hModel, *(varNames+[selND.wName]))
        else:
            plotDF = plotFun(hModel, *varNames)
        self.plotResults[plot.name] = plotDF

    @staticmethod
    def makePlotModel(plot):
        import ROOT
        modCls = getattr(ROOT.RDF, "TH{0:d}DModel".format(len(plot.binnings)))
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
