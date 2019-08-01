"""
ROOT::RDataFrame backend classes
"""
import logging
logger = logging.getLogger(__name__)

import copy
from itertools import chain, count
from functools import partial

from .plots import FactoryBackend, Selection
from . import treefunctions as op
from . import treeoperations as top

class SelWithDefines(top.CppStrRedir):
    def __init__(self, parent, variation="nominal"):
        self.explDefine = list()
        self.var = None ## nodes for related systematic variations (if different by more than the weight)
        if isinstance(parent, SelWithDefines):
            self.parent = parent
            self.backend = parent.backend
            self.df = parent.df
            self._definedColumns = dict(parent._definedColumns)
            if variation == "nominal":
                self.var = dict((varNm, SelWithDefines(pvar, pvar.df, variation=varNm)) for varNm, pvar in parent.var.items())
                ## parent is either nominal (when branching off), or the corresponding variation of the parent of the nominal node
        elif isinstance(parent, DataframeBackend):
            self.parent = None
            self.backend = parent
            self.df = parent.rootDF
            self._definedColumns = dict()
            self.var = dict()
        else:
            raise RuntimeError("Can only define SelWithDefines with a DataframeBackend or another SelWithDefines")
        self.wName = dict()
        top.CppStrRedir.__init__(self)

    def addWeight(self, weights=None, wName=None, variation="nominal"):
        parentWeight = None
        if self.parent and self.parent.wName:
            parentWeight = self.parent.wName.get(variation, self.parent.wName["nominal"])
        if weights:
            weightExpr = Selection._makeExprProduct(
                ([top.adaptArg(op.extVar("float", parentWeight), typeHint="float")]+weights) if parentWeight
                else weights
                )
            self._define(wName, weightExpr)
            self.wName[variation] = wName
        elif parentWeight:
            self.wName[variation] = parentWeight
        else:
            assert not wName
            self.wName[variation] = None

    def addCut(self, cuts):
        cutExpr = Selection._makeExprAnd(cuts)
        cutStr = cutExpr.get_cppStr(defCache=self)
        self._addFilterStr(cutStr)

    def _addFilterStr(self, filterStr): ## add filter with string already made
        logger.debug("Filtering with {0}".format(filterStr))
        self.df = self.df.Filter(filterStr)

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
        self.selDFs = dict()      ## (selection name, variation) -> SelWithDefines
        self.plotResults = dict() ## plot name -> list of result pointers
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
        cutStr = None
        if sele._cuts:
            assert sele.parent ## FIXME there *needs* to be a root no-op sel
            ## trick: by passing defCache=parentDF and doing this *before* constructing the nominal node,
            ## any definitions end up in the node above, and are in principle available for other sub-selections too
            cutStr = Selection._makeExprAnd(sele._cuts).get_cppStr(defCache=self.selDFs[sele.parent.name])
        selnd = SelWithDefines(self.selDFs[sele.parent.name] if sele.parent else self)
        if cutStr:
            selnd._addFilterStr(cutStr)
        selnd.addWeight(weights=sele._weights, wName=("w_{0}".format(sele.name) if sele._weights else None))
        ## Next: loop through all systematics that affect the cuts or weights
        ## - variation is there (i.e. it affected previous cuts) -> add cut and weight to var
        ## - variation is not there:
        ##   - only weight -> keep going as before
        ##   - cut -> create variations
        ## - with systName and variations, the distinction between collection and scalefactor is irrelevant
        ## -> make them derive from the same base class that defines this interface

        for systN,systVars in sele.weightSystematics.items(): ## weights-only systematics (do not need a new node, nominal selection)
            logger.debug("Adding weight variations {0} for systematic {1}".format(systVars, systN))
            wfToChange = []
            wfKeep = list(sele._weights)
            isthissyst = partial((lambda sN,iw : isinstance(iw, top.ScaleFactorWithSystOp) and iw.systName == sN), systN)
            if systN in sele._wSysts:
                for wf in sele._weights:
                    if any(top.collectNodes(wf, select=isthissyst)):
                        wfToChange.append(wf)
                        del wfKeep[wfKeep.index(wf)]
            for vard in systVars:
                varn = "{0}{1}".format(systN, vard)
                wfChanged = []
                for wf in wfToChange:
                    newf = copy.deepcopy(wf)
                    for nd in top.collectNodes(newf, select=isthissyst):
                        nd.changeVariation(vard)
                    wfChanged.append(newf)
                selnd.addWeight(weights=(wfKeep+wfChanged), wName=("w_{0}__{1}".format(sele.name, varn) if sele._weights else None), variation=varn)
        self.selDFs[sele.name] = selnd
        ## TODO add non-weightonly systematics (will need new nodes; will be: cuts and/or weights depend on collections)

    def addPlot(self, plot):
        """ Define ROOT::RDataFrame objects needed for this plot (and keep track of the result pointer) """
        if plot.name in self.plotResults:
            raise ValueError("A Plot with the name '{0}' already exists".format(plot.name))
        varSysts = dict((sfs.systName, sfs.variations) for sfs in chain.from_iterable(
            top.collectNodes(vi, select=(lambda nd : isinstance(nd, top.SystModifiedCollectionOp) and nd.systName and len(nd.variations) > 1))
            for vi in plot.variables))
        if varSysts:
            logger.debug("Plot variables are affected by systematics {0!s}".format(varSysts))
        ## TODO DataFrame might throw (or segfault), but we should catch all possible errors before
        ## NOTE pre/postfixes should already be inside the plot name
        selND = self.selDFs[plot.selection.name] ## TODO update with variations (later)
        varExprs = dict(("v{0:d}_{1}".format(i, plot.name), selND(var)) for i,var in zip(count(), plot.variables))
        plotDF = selND.df ## after getting the expressions, to pick up columns that were defined on-demand
        for vName, vExpr in varExprs.items():
            #logger.debug("Defining {0} as {1} (defined for {2}: {3})".format(vName, vExpr, plotDF, ", ".join(plotDF.GetDefinedColumnNames()))) ## needs 6.16
            logger.debug("Defining {0} as {1}".format(vName, vExpr))
            plotDF = plotDF.Define(vName, vExpr)
        plotFun = getattr(plotDF, "Histo{0:d}D".format(len(plot.variables)))
        ## TODO for weight-only: essentially fill for all weights (incl. nominal)
        ## TODO: the other cases will be
        ## - cut depends but not plotvar -> repeat on different selnd with nomianl
        ## - cut does not depends but plotvar -> repeat on the same selnd
        ## - both cut and plotvar depend -> run for matching combinations
        ## (in all cases, weight *may* depend or not)
        if selND.wName["nominal"]: ## nontrivial weight
            logger.debug("Adding plot {0} with variables {1} and weights {2}".format(plot.name, ", ".join(varExprs.keys()), ", ".join(selND.wName.values())))
            plotDF = [ plotFun( DataframeBackend.makePlotModel(plot, variation=wvarName),
                                *chain(varExprs.keys(), [ varWeight ]) )
                        for wvarName, varWeight in selND.wName.items() ]
        else:
            logger.debug("Adding plot {0} with variables {1}".format(plot.name, ", ".join(varExprs.keys())))
            hModel = DataframeBackend.makePlotModel(plot)
            plotDF = [ plotFun(hModel, *varExprs.keys()) ]
        self.plotResults[plot.name] = plotDF

    @staticmethod
    def makePlotModel(plot, variation="nominal"):
        from cppyy import gbl
        modCls = getattr(gbl.RDF, "TH{0:d}DModel".format(len(plot.binnings)))
        name = plot.name
        if variation != "nominal":
            name = "__".join((name, variation))
        return modCls(name, plot.title, *chain.from_iterable(
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

    def getPlotResults(self, plot):
        return self.plotResults[plot.name]

    def writeSkim(self, sele, outputFile, treeName, definedBranches=None, origBranchesToKeep=None, maxSelected=-1):
        selND = self.selDFs[sele.name]

        allcolN = selND.df.GetColumnNames()
        defcolN = selND.df.GetDefinedColumnNames()
        colNToKeep = type(allcolN)()
        if origBranchesToKeep is None: ## keep all if not defined
            for cn in allcolN:
                if cn not in defcolN:
                    colNToKeep.push_back(cn)
        elif len(origBranchesToKeep) != 0:
            for cn in origBranchesToKeep:
                if cn not in allcolN:
                    raise RuntimeError("Requested column '{0}' from input not found".format(cn))
                if cn in defcolN:
                    raise RuntimeError("Requested column '{0}' from input is a defined column".format(cn))
                colNToKeep.push_back(cn)

        for dN, dExpr in definedBranches.items():
            selND._define(dN, top.adaptArg(dExpr))
            colNToKeep.push_back(dN)

        selDF = selND.df
        if maxSelected > 0:
            selDF = selDF.Range(maxSelected)

        selDF.Snapshot(treeName, outputFile, colNToKeep)
