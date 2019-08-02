"""
ROOT::RDataFrame backend classes
"""
import logging
logger = logging.getLogger(__name__)

import copy
from itertools import chain
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

    def addWeight(self, weights=None, wName=None, parentWeight=None, variation="nominal"):
        if variation == "nominal" and parentWeight is None and self.parent and self.parent.wName:
            parentWeight = self.parent.wName["nominal"]
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

        super(SystModifiedCollectionOp, self).__init__(wrapped, name)
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
        nomParentNd = self.selDFs[sele.parent.name] if sele.parent else None
        if sele._cuts:
            assert sele.parent ## FIXME there *needs* to be a root no-op sel
            ## trick: by passing defCache=parentDF and doing this *before* constructing the nominal node,
            ## any definitions end up in the node above, and are in principle available for other sub-selections too
            cutStr = Selection._makeExprAnd(sele._cuts).get_cppStr(defCache=self.selDFs[sele.parent.name])
        nomNd = SelWithDefines(nomParentNd if nomParentNd else self)
        if cutStr:
            nomNd._addFilterStr(cutStr)
        nomNd.addWeight(weights=sele._weights, wName=("w_{0}".format(sele.name) if sele._weights else None))
        self.selDFs[sele.name] = nomNd

        ## Next: loop through all systematics that affect the cuts or weights
        ## - variation is there (i.e. it affected previous cuts) -> add cut and weight to var
        ## - variation is not there:
        ##   - only weight -> keep going as before
        ##   - cut -> create variations
        ## - with systName and variations, the distinction between collection and scalefactor is irrelevant
        ## -> make them derive from the same base class that defines this interface

        weightSyst = sele.weightSystematics
        cutSyst = sele.cutSystematics
        for systN, systVars in sele.systematics.items(): ## the two above merged
            logger.debug("Adding weight variations {0} for systematic {1}".format(systVars, systN))
            ## figure out which cuts and weight factors are affected by this systematic
            isthissyst = partial((lambda sN,iw : isinstance(iw, top.OpWithSyst) and iw.systName == sN), systN)
            ctToChange = []
            ctKeep = list(sele._cuts)
            wfToChange = []
            wfKeep = list(sele._weights)
            if systN in cutSyst:
                nRem = 0
                for i,ct in enumerate(sele._cuts):
                    if any(top.collectNodes(ct, select=isthissyst)):
                        ctToChange.append(ct)
                        del ctKeep[i-nRem]
                        nRem += 1
            if systN in weightSyst:
                nRem = 0
                for i,wf in enumerate(sele._weights):
                    if any(top.collectNodes(wf, select=isthissyst)):
                        wfToChange.append(wf)
                        del wfKeep[i-nRem]
                        nRem += 1
            ## construct variation nodes (if necessary)
            for vard in systVars:
                varn = "{0}{1}".format(systN, vard)
                ## add cuts to the appropriate node, if affected by systematics (here or up)
                varParentNd = None ## set parent node if not the nominal one
                if nomParentNd and varn in nomParentNd.var: ## -> continue on branch
                    varParentNd = nomParentNd.var[varn]
                elif ctToChange: ## -> branch off now
                    varParentNd = nomParentNd
                if not varParentNd: ## cuts unaffected (here and in parent), can stick with nominal
                    varNd = nomNd
                else: ## on branch, so add cuts (if any)
                    if len(sele._cuts) == 0: ## no cuts, reuse parent
                        varNd = varParentNd
                    else:
                        ctChanged = []
                        for ct in ctToChange: ## empty if sele._cuts are not affected
                            newct = copy.deepcopy(ct)
                            for nd in top.collectNodes(newct, select=isthissyst):
                                nd.changeVariation(vard)
                            ctChanged.append(newct)
                        cutStr = Selection._makeExprAnd(ctKeep+ctChanged).get_cppStr(defCache=varParentNd)
                        varNd = SelWithDefines(varParentNd)
                        varNd._addFilterStr(cutStr)
                    nomNd.var[varn] = varNd
                ## next: attach weights (modified if needed) to varNd
                if varParentNd:
                    parwn = varParentNd.wName.get(varn, varParentNd.wName.get("nominal"))
                elif nomParentNd:
                    parwn = nomParentNd.wName.get(varn, nomParentNd.wName.get("nominal"))
                else:
                    parwn = None ## no prior weights at all
                if not sele._weights:
                    logger.debug("{0} systematic variation {1}: reusing {2}".format(sele.name, varn, parwn))
                    varNd.addWeight(parentWeight=parwn, variation=varn)
                else:
                    if wfToChange or varNd != nomNd or ( nomParentNd and varn in nomParentNd.wName ):
                        wfChanged = []
                        for wf in wfToChange:
                            newf = copy.deepcopy(wf)
                            for nd in top.collectNodes(newf, select=isthissyst):
                                nd.changeVariation(vard)
                            wfChanged.append(newf)
                        logger.debug("{0} systematic variation {1}: defining new weight based on {2}".format(sele.name, varn, parwn))
                        varNd.addWeight(weights=(wfKeep+wfChanged), wName=("w_{0}__{1}".format(sele.name, varn) if sele._weights else None), parentWeight=parwn, variation=varn)
                    else: ## varNd == nomNd, not branched, and parent does not have weight variation
                        logger.debug("{0} systematic variation {1}: reusing nominal {2}".format(sele.name, varn, varNd.wName["nominal"]))
                        varNd.addWeight(parentWeight=varNd.wName["nominal"], variation=varn)

    def addPlot(self, plot):
        """ Define ROOT::RDataFrame objects needed for this plot (and keep track of the result pointer) """
        if plot.name in self.plotResults:
            raise ValueError("A Plot with the name '{0}' already exists".format(plot.name))
        varSysts = dict((sfs.systName, sfs.variations) for sfs in chain.from_iterable(
            top.collectNodes(vi, select=(lambda nd : isinstance(nd, top.OpWithSyst) and nd.systName and nd.variations))
            for vi in plot.variables))
        if varSysts:
            logger.debug("Plot variables are affected by systematics {0!s}".format(varSysts))
        selSysts = plot.selection.systematics
        allSysts = dict(plot.selection.systematics)
        allSysts.update(varSysts)
        if allSysts:
            logger.debug("Plot is affected by systematics {0!s} through selection or variable(s)".format(allSysts))

        nomNd = self.selDFs[plot.selection.name]
        plotRes = []
        ## Add nominal plot
        nomVarExprs = dict(("v{0:d}_{1}".format(i, plot.name), nomNd(var)) for i,var in enumerate(plot.variables))
        nomPlotDF = nomNd.df
        for vName, vExpr in nomVarExprs.items():
            logger.debug("Defining {0} as {1}".format(vName, vExpr))
            nomPlotDF = nomPlotDF.Define(vName, vExpr)
        nomPlotFun = getattr(nomPlotDF, "Histo{0:d}D".format(len(plot.variables)))
        plotModel = DataframeBackend.makePlotModel(plot)
        if nomNd.wName["nominal"]: ## nontrivial weight
            logger.debug("Adding plot {0} with variables {1} and weight {2}".format(plot.name, ", ".join(nomVarExprs.keys()), nomNd.wName["nominal"]))
            plotRes.append(nomPlotFun(plotModel, *chain(nomVarExprs.keys(), [ nomNd.wName["nominal"] ]) ))
        else: # no weight
            logger.debug("Adding plot {0} with variables {1}".format(plot.name, ", ".join(nomVarExprs.keys())))
            plotRes.append(nomPlotFun(plotModel, *nomVarExprs.keys()))

        ## Same for all the systematics
        for systN, systVars in allSysts.items():
            isthissyst = partial((lambda sN,iw : isinstance(iw, top.OpWithSyst) and iw.systName == sN), systN)
            idxVarsToChange = []
            for i,xvar in enumerate(plot.variables):
                if any(top.collectNodes(xvar, select=isthissyst)):
                    idxVarsToChange.append(i)
            for vard in systVars:
                varn = "{0}{1}".format(systN, vard)
                if systN in varSysts or varn in nomNd.var:
                    varNd = nomNd.var.get(varn, nomNd)
                    varExprs = {}
                    for i,xvar in enumerate(plot.variables):
                        if i in idxVarsToChange:
                            varVar = copy.deepcopy(xvar)
                            for nd in top.collectNodes(varVar, select=isthissyst):
                                nd.changeVariation(vard)
                        else:
                            varVar = xvar
                        varExprs["v{0:d}_{1}__{2}".format(i, plot.name, varn)] = varNd(varVar)
                    plotDF = varNd.df
                    for vName, vExpr in varExprs.items():
                        logger.debug("Defining {0} as {1}".format(vName, vExpr))
                        plotDF = plotDF.Define(vName, vExpr)
                    plotFun = getattr(plotDF, "Histo{0:d}D".format(len(plot.variables)))
                    plotModel = DataframeBackend.makePlotModel(plot, variation=varn)
                    if varn not in varNd.wName:
                        logger.error("{0} not in {1!s}".format(varn, varNd.wName.keys()))
                    wN = varNd.wName[varn] if systN in selSysts else varNd.wName["nominal"] ## else should be "only in the variables", so varNd == nomNd then
                    if wN is not None: ## nontrivial weight
                        logger.debug("Adding variation {0} plot {1} with variables {2} and weight {3}".format(varn, plot.name, ", ".join(varExprs.keys()), wN))
                        plotRes.append(plotFun(plotModel, *chain(varExprs.keys(), [ wN ]) ))
                    else: ## no weight
                        logger.debug("Adding variation {0} plot {1} with variables {2}".format(varn, plot.name, ", ".join(varExprs.keys())))
                        plotRes.append(plotFun(plotModel, *varExprs.keys()))
                else: ## can reuse variables, but may need to take care of weight
                    wN = nomNd.wName[varn]
                    if wN is not None:
                        logger.debug("Adding variation {0} plot {1} with variables {2} and weight {3}".format(varn, plot.name, ", ".join(nomVarExprs.keys()), wN))
                        plotRes.append(nomPlotFun(plotModel, *chain(nomVarExprs.keys(), [ wN ]) ))
                    else: ## no weight
                        logger.error("A systematic that doesn't affect cuts, variables, or a weight... this is weird ")
                        logger.debug("Adding variation {0} plot {1} with variables {2}".format(varn, plot.name, ", ".join(nomVarExprs.keys())))
                        plotRes.append(nomPlotFun(plotModel, *nomVarExprs.keys()))
        ## at the end
        self.plotResults[plot.name] = plotRes

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
