"""
Main user-visible module of tree operations, separate from the operation definitions
and proxy classes
User-facing module
"""

from collections import defaultdict
from functools import partial

from .treeproxies import *
from . import treefunctions as op

import logging
logger = logging.getLogger(__name__)

def allLeafs(branch):
    """
    Recursively collect TTree leaves (TLeaf and TBranchElement)
    """
    for br in branch.GetListOfBranches():
        from .root import gbl
        if isinstance(br, gbl.TBranchElement):
            yield br
        else:
            for lv in br.GetListOfLeaves():
                yield lv

def normVarName(varName):
    """ Normalize variation name: if ending in up or down, make sure this part has no capitals (for plotIt) """
    if len(varName) >= 2 and varName[-2:].upper() == "UP":
        return "{0}up".format(varName[:-2])
    elif len(varName) >= 4 and varName[-4:].upper() == "DOWN":
        return "{0}down".format(varName[:-4])
    else:
        return varName

def getMETVarName_calc(nm):
    if nm in ("pt", "phi"):
        return nm, "raw"

def getJetEnergyVarName_calc(nm):
    if nm in ("pt", "mass"):
        return nm, "raw"
def getMuonMomentumVarName_calc(nm):
    if nm == "pt":
        return nm, "raw"

## TODO: make configurable
def getJetMETVarName_postproc(nm, collgrpname=None):
    if nm.split("_")[0] in ("pt", "eta", "phi", "mass") and len(nm.split("_")) >= 2:
        return (nm.split("_")[0], "_".join(nm.split("_")[1:]))
    elif nm.startswith("btagSF"):
        for tagger in ["csvv2", "deepcsv", "deepjet", "cmva"]:
            for wp in ["L", "M", "T", "shape"]:
                sfName = f"btagSF_{tagger}_{wp}"
                if not nm.startswith(sfName):
                    continue
                if nm == sfName:
                    return nm, "nom"
                upOrDown = "up" if "up" in nm else "down"
                if wp != "shape": # b-tag SF systematics
                    return sfName, f"{sfName}{upOrDown}"
                else:
                    syst = nm.split(f"{sfName}_{upOrDown}_")[1]
                    if "jes" not in syst: # b-tag shape systematics
                        return sfName, f"{sfName}_{syst}{upOrDown}"
                    else: # jes systematics
                        if syst == "jes":
                            syst = "jesTotal"
                        else:
                            return sfName, f"{syst}{upOrDown}"

## Attribute classes (like property) to customize the proxy classes

class proxy(object): ## the default one
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op.result

class funProxy(object): ## the generic one
    def __init__(self, fun):
        self.fun = fun
    def __get__(self, inst, cls):
        return self.fun(inst)

class itemProxy(object):
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op[inst._idx]

class itemRefProxy(object):
    def __init__(self, op, getTarget):
        self.op = op
        self.getTarget = getTarget
    def __get__(self, inst, cls):
        return self.getTarget(inst)[self.op[inst._idx]]

class itemObjProxy(object): ## re-construct an object that was split in arrays
    def __init__(self, typeName, args):
        self.typeName = typeName
        self.args = args
    def __get__(self, inst, cls):
        return Construct(self.typeName, tuple(arg[inst._idx] for arg in self.args)).result

class altProxy(object): ## grouped: instance has the branch map
    def __init__(self, name, op):
        self.name = name
        self.op = op
    def __get__(self, inst, cls):
        return inst.brMap.get(self.name, self.op).result

class altItemProxy(object): ## collection carries branch map (instance comes from collection.__getitem__)
    def __init__(self, name, op):
        self.name = name
        self.op = op
    def __get__(self, inst, cls):
        return inst._parent.brMap.get(self.name, self.op).result[inst._idx]

def decorateTTW(aTree, description=None):
    ## NOTE: WORK IN PROGRESS
    if description is None:
        description = dict()

    class SetAsParentOfAtt(object): ## FIXME can be removed in favour of below, probably
        ## callback to set the tree as _parent of the groups
        def __init__(self, name):
            self.name = name
        def __call__(self, inst):
            getattr(inst, self.name)._parent = inst

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format("ttW")}
    ## fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm))) for lvNm,lv in allTreeLeafs.items()))
    tree_postconstr = []
    for grpNm, desc in description.items():
        if desc["type"] == "container":
            prefix = desc["prefix"]
            desc_rem = dict(desc)
            desc_rem.pop("type")
            desc_rem.pop("prefix")
            itm_dict = {
                "__doc__" : "{0} proxy class".format(grpNm)
                }
            itm_lvs_vec = set()
            itm_lvs_arrCnt = dict()
            from .root import gbl
            for lvNm, lv in allTreeLeafs.items():
                if lvNm.startswith(prefix):
                    if isinstance(lv, gbl.TBranchElement):
                        itm_lvs_vec.add(lvNm)
                    else:
                        if lv.GetLeafCount():
                            cntLvNm = lv.GetLeafCount().GetName()
                            if cntLvNm not in itm_lvs_arrCnt:
                                itm_lvs_arrCnt[lvNm] = []
                            itm_lvs_arrCnt[cntLvNm].append(lvNm)
                        else:
                            if lvNm not in itm_lvs_arrCnt:
                                itm_lvs_arrCnt[lvNm] = []
            for cntLvNm,arrLvs in itm_lvs_arrCnt.items():
                sizeOp = GetColumn(allTreeLeafs[cntLvNm].GetTypeName(), cntLvNm)
                if cntLvNm.endswith("_") and cntLvNm.rstrip("_") in itm_lvs_vec:
                    if all(aLvNm[len(cntLvNm):].startswith("fCoordinates.f") for aLvNm in arrLvs) and len(arrLvs) == 4:
                        logger.debug("Detected a split LorentzVector array {}".format(cntLvNm.rstrip("_")))
                        atts = [ aLvNm[len(cntLvNm)+len("fCoordinates.f"):] for aLvNm in arrLvs ]
                        logger.debug(" ".join(aLvNm[len(cntLvNm)+len("fCoordinates.f"):] for aLvNm in arrLvs))
                        tpNm = next(allTreeLeafs[aLvNm].GetTypeName() for aLvNm in arrLvs)
                        itm_dict[cntLvNm[len(prefix):].rstrip("_")] = itemObjProxy("ROOT::Math::LorentzVector<ROOT::Math::{0}4D<{1}> >".format("".join(atts), tpNm),
                                tuple(GetArrayColumn(allTreeLeafs[aLvNm].GetTypeName(), aLvNm, sizeOp).result for aLvNm in arrLvs))
                        for aLvNm in arrLvs:
                            del tree_dict[aLvNm]
                        del tree_dict[cntLvNm]
                        del tree_dict[cntLvNm.rstrip("_")]
                        itm_lvs_vec.remove(cntLvNm.rstrip("_"))
                    elif all(aLvNm[len(cntLvNm):].startswith("fCoordinates.f") for aLvNm in arrLvs) and len(arrLvs) == 3:
                        atts = [ aLvNm[len(cntLvNm)+len("fCoordinates.f"):] for aLvNm in arrLvs ]
                        assert("X" in atts and "Y" in atts and "Z" in atts)
                        tpNm = next(allTreeLeafs[aLvNm].GetTypeName() for aLvNm in arrLvs)
                        itm_dict[cntLvNm[len(prefix):].rstrip("_")] = itemObjProxy("ROOT::Math::PositionVector3D<ROOT::Math::Cartesian3D<{0}>,ROOT::Math::DefaultCoordinateSystemTag>".format(tpNm),
                                tuple(GetArrayColumn(allTreeLeafs[aLvNm].GetTypeName(), aLvNm, sizeOp).result for alvNm in arrLvs))
                        for aLvNm in arrLvs:
                            del tree_dict[aLvNm]
                        del tree_dict[cntLvNm]
                        del tree_dict[cntLvNm.rstrip("_")]
                        itm_lvs_vec.remove(cntLvNm.rstrip("_"))
                    else:
                        logger.warning("Detected a split object that's not a LorentzVector. Attributes are: {}".format(", ".join(aLvNm[len(cntLvNm):] for aLvNm in arrLvs)))
                else:
                    for arrLv in arrLvs:
                        lvNm_short = arrLv[len(prefix):]
                        col = GetArrayColumn(allTreeLeafs[arrLv].GetTypeName(), arrLv, sizeOp).result
                        itm_dict[lvNm_short] = itemProxy(col)
                        del tree_dict[arrLv]
                    del tree_dict[cntLvNm]
            aLvNm = next(lv for lv in itm_lvs_vec if "{0}_".format(lv) not in itm_lvs_vec)
            sizeOp = adaptArg(op.rng_len(GetColumn(allTreeLeafs[aLvNm].GetClassName(), aLvNm).result)) ## abuse a bit
            for lvNm in itm_lvs_vec:
                lvNm_short = lvNm[len(prefix):]
                col = GetColumn(allTreeLeafs[lvNm].GetClassName(), lvNm).result
                itm_dict[lvNm_short] = itemProxy(col)
                del tree_dict[lvNm]
            itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
            tree_dict[grpNm] = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            tree_postconstr.append(SetAsParentOfAtt(grpNm))
            logger.info("Processed group {0}, remaining description: {1}".format(grpNm, desc_rem))
        elif desc["type"] == "group":
            prefix = desc["prefix"]
            desc_rem = dict(desc)
            desc_rem.pop("type")
            desc_rem.pop("prefix")
            grp_dict = {
                "__doc__" : "{0} leaf group proxy class".format(grpNm)
                }
            grp_lvNms = set(lvNm for lvNm in allTreeLeafs.keys() if lvNm.startswith(prefix))
            grp_dict.update(dict((lvNm[len(prefix):], proxy(GetColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm))) for lvNm in grp_lvNms))
            grpcls = type("{0}LeafGroupProxy".format(grpNm), (LeafGroupProxy,), grp_dict)
            tree_dict[grpNm] = grpcls(grpNm, None) ## NOTE set later
            for lvNm in grp_lvNms:
                del tree_dict[lvNm]
            tree_postconstr.append(SetAsParentOfAtt(grpNm))

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)

    for pc in tree_postconstr:
        pc(treeProxy)

    return treeProxy

def _makeAltClassAndMaps(name, dict_orig, getVarName, systName=None, nomName="nom", exclVars=None,
        getCol=lambda op : op, attCls=None, altBases=None): ## internal, leaf/group/collection
    ## getVarName should return the variable and variation name (nomName for the nominal one)
    ## if this is a systematic variation branch - otherwise None
    dict_alt = dict(dict_orig)
    ## collect ops of kinematic variables that change (nominal as well as varied)
    var_atts = defaultdict(dict)
    for nm, nmAtt in dict_orig.items():
        test = getVarName(nm, collgrpname=name)
        if test is not None:
            attNm, varNm = test
            var_atts[attNm][normVarName(varNm)] = nmAtt.op
            del dict_alt[nm]
    ## redirect in altProxy
    for attNm in var_atts.keys():
        dict_alt[attNm] = attCls(attNm, dict_orig[attNm].op)
    cls_alt = type("Alt{0}Proxy".format(name), altBases, dict_alt)
    ## construct the map of maps of redirections { variation : { varName : op } }
    brMapMap = {}
    for attNm,vAtts in var_atts.items():
        for var,vop in vAtts.items():
            if var not in brMapMap:
                brMapMap[var] = {}
            brMapMap[var][attNm] = getCol(vop)
    ## nominal: with systematic variations (all are valid, but not all need to modify)
    allVars = list(k for k in brMapMap.keys() if k not in exclVars and k != nomName)
    brMapMap["nomWithSyst"] = dict((attNm,
        SystAltOp(
            getCol(vAtts[nomName]), systName,
            dict((var, getCol(vop)) for var,vop in vAtts.items() if var not in exclVars),
            valid=tuple(var for var in allVars if var in vAtts),
            ))
        for attNm,vAtts in var_atts.items())
    return cls_alt, brMapMap

## Helper classes
class NanoSystematicVarSpec:
    """ Interface for classes that specify how to incorporate systematics or on-the-fly corrections in the decorated tree """
    def __init__(self, systName, nomName="nom", origName=None, exclVars=None, isCalc=False):
        self.systName = systName
        self.nomName = nomName
        self.origName = origName
        self.exclVars = exclVars if exclVars is not None else tuple()
        self.isCalc = isCalc
    def appliesTo(self, name):
        """ Return true if this systematic variation requires action for this variable, group, or collection """
        return False
    def getVarName(self, branchName, collgrpname=None):
        pass

class ReadVariableVarWithSuffix(NanoSystematicVarSpec):
    """ Read variations of a single branch from branches with the same name with a suffix """
    def __init__(self, commonName, sep="_", systName=None, nomName="nom", exclVars=None):
        super(ReadVariableVarWithSuffix, self).__init__(
                (systName if systName is not None else commonName),
                nomName=nomName, exclVars=exclVars, isCalc=False)
        self.prefix = commonName
        self.sep = sep
    def appliesTo(self, name):
        return name.startswith(self.prefix)
    def getVarName(self, branchName, collgrpname=None):
        variNm = normVarName(branchName[len(self.prefix):].lstrip(self.sep))
        return self.prefix, variNm if variNm else self.nomName

nanoPUWeightVar = ReadVariableVarWithSuffix("puWeight")

class ReadJetMETVar(NanoSystematicVarSpec):
    def __init__(self, jetsName, jetAttrs, metName, metAttrs, systName="jet", nomName="nom", origName="raw", exclVars=None):
        super(ReadJetMETVar, self).__init__(
                systName, nomName=nomName, origName=origName, isCalc=True,
                exclVars=(exclVars if exclVars is not None else (origName,)))
        self.jetsName = jetsName
        self.metName = metName
        self.jetAttrs = jetAttrs
        self.metAttrs = metAttrs
    def appliesTo(self, name):
        return name in (self.jetsName, self.metName)
    def getVarName(self, nm, collgrpname=None):
        pass ## TODO migrate above here

class CalcCollection(NanoSystematicVarSpec):
    def __init__(self, collName, calcAttrs, systName=None, nomName="nominal", origName="raw", exclVars=None):
        super(CalcCollection, self).__init__(
                (systName if systName is not None else collName.lower()),
                exclVars=(exclVars if exclVars is not None else (origName,)),
                nomName=nomName, origName=origName, isCalc=True)
        self.collName = collName
        self.calcAttrs = calcAttrs
    def appliesTo(self, name):
        return name == self.collName
    def getVarName(self, nm, collgrpname=None):
        if nm in self.calcAttrs:
            return nm, self.origName

class CalcJetMETVar(NanoSystematicVarSpec):
    def __init__(self, jetsName, jetAttrs, metName, metAttrs, systName="jet", nomName="nominal", origName="raw", exclVars=None):
        super(CalcJetMETVar, self).__init__(
                systName, nomName=nomName, origName=origName, isCalc=True,
                exclVars=(exclVars if exclVars is not None else (origName,)))
        self.jetsName = jetsName
        self.metName = metName
        self.jetAttrs = jetAttrs
        self.metAttrs = metAttrs
    def appliesTo(self, name):
        return name in (self.jetsName, self.metName)
    def getVarName(self, nm, collgrpname=None):
        if collgrpname == self.jetsName and nm in self.jetAttrs:
            return nm, self.origName
        if collgrpname == self.metName and nm in self.metAttrs:
            return nm, self.origName

nanoRochesterCalc = CalcCollection("Muon", ("pt",))
nanoJetMETCalc = CalcJetMETVar("Jet", ("pt", "mass"), "MET", ("pt", "phi"))
nanoJetMETCalc_METFixEE2017 = CalcJetMETVar("Jet", ("pt", "mass"), "METFixEE2017", ("pt", "phi"))

def decorateNanoAOD(aTree, description=None, isMC=False, systVariations=[ nanoPUWeightVar ]):
    """ Decorate a CMS NanoAOD Events tree

    Variation branches following the NanoAODTools conventions (e.g. Jet_pt_nom)
    are automatically used (but calculators for the same collection take
    precendence, if requested).

    :param aTree: TTree to decorate
    :param description: reserved for future version/flavour-dependence
    :param isMC: simulation or not
    :param systVariations: systematic variations and on-th-fly corrections to add (list of :py:class:`bamboo.treedecorators.NanoSystematicVarSpec` instances)
    """
    if description is None:
        description = dict()

    class GetItemRefCollection:
        def __init__(self, name):
            self.name = name
            self._parent = None
        def __call__(self, me):
            return getattr(self._parent, self.name)
    class GetItemRefCollection_toVar:
        def __init__(self, name):
            self.name = name
            self._parent = None
        def __call__(self, me):
            return getattr(self._parent, self.name).orig

    def addP4ToObj(prefix, lvNms):
        return funProxy(partial( (lambda getEta,getM,inst:
                   Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >",
                             (inst.pt, getEta(inst), inst.phi, getM(inst))).result),
                   ((lambda inst: inst.eta) if f"{prefix}eta" in lvNms else (lambda inst: 0.)),
                   ((lambda inst: inst.mass) if f"{prefix}mass" in lvNms else (lambda inst: 0.)) ))

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    ## NOTE first attempt: fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm) if not lv.GetLeafCount()
        else GetArrayColumn(lv.GetTypeName(), lvNm, GetColumn(allTreeLeafs[lv.GetLeafCount().GetName()].GetTypeName(), lv.GetLeafCount().GetName()))
        )) for lvNm,lv in allTreeLeafs.items()))
    tree_children = list()
    def setTreeAtt(name, proxy, setParent=True):
        tree_dict[name] = proxy
        if setParent:
            tree_children.append(proxy)

    ## variables with variations
    for vari in systVariations:
        toRem = []
        brMap = {}
        attNm = None
        for nm,nmAtt in tree_dict.items():
            if vari.appliesTo(nm): ## only for single-variable variations
                attNm,varNm = vari.getVarName(nm)
                brMap[varNm] = nmAtt.op
                toRem.append(nm)
        if brMap and ( len(brMap) > 1 or vari.isCalc ):
            for nm in toRem:
                del tree_dict[nm]
            logger.debug(f"Detected systematic variations for variable {attNm} (variations: {list(brMap.keys())!r}")
            brMap["nomWithSyst"] = SystAltOp(brMap[vari.nomName], vari.systName,
                dict(brMap), valid=tuple(var for var in brMap.keys() if var != vari.nomName)
                )
            varsProxy = AltLeafVariations(None, brMap, typeName=brMap[vari.nomName].typeName)
            setTreeAtt(f"_{attNm}", varsProxy)
            setTreeAtt(attNm, varsProxy["nomWithSyst"], False)

    ## non-collection branches to group
    simpleGroupPrefixes = ("CaloMET_", "ChsMET_", "MET_", "PV_", "PuppiMET_", "RawMET_", "TkMET_", "Flag_", "HLT_", "L1_") ## TODO get this from description?
    simpleGroupPrefixes_opt = ("METFixEE2017_",)
    simpleGroupPrefixes_Gen = ("GenMET_", "Generator_", "LHE_", "HTXS_")
    ## check which are there, and for which we need to read variations
    grp_found = []
    grp_readVar = []
    for prefix in (chain(simpleGroupPrefixes, simpleGroupPrefixes_opt, simpleGroupPrefixes_Gen) if isMC else chain(simpleGroupPrefixes, simpleGroupPrefixes_opt)):
        if prefix not in simpleGroupPrefixes_opt and not any(lvNm.startswith(prefix) for lvNm in allTreeLeafs):
            logger.warning("No branch name starting with {0} in the tree - skipping group".format(prefix))
        else:
            grp_found.append(prefix)
            if prefix.startswith("MET") and "{0}_pt_nom".format(prefix.rstrip("_")) in allTreeLeafs:
                grp_readVar.append(prefix)
    for prefix in grp_found:
        grpNm = prefix.rstrip("_")
        grp_dict = {
            "__doc__" : "{0} leaf group proxy class".format(grpNm)
            }
        grp_lvNms = set(lvNm for lvNm in allTreeLeafs.keys() if lvNm.startswith(prefix))
        grp_dict.update(dict((lvNm[len(prefix):], proxy(GetColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm))) for lvNm in grp_lvNms))
        if f"{prefix}pt" in grp_lvNms and f"{prefix}phi" in grp_lvNms:
            grp_dict["p4"] = addP4ToObj(prefix, grp_lvNms)
        grpcls = type("{0}LeafGroupProxy".format(grpNm), (LeafGroupProxy,), grp_dict)
        for lvNm in grp_lvNms:
            del tree_dict[lvNm]
        ## default group proxy, replaced below if needed
        grp_proxy = grpcls(grpNm, None)
        setTreeAtt(grpNm, grp_proxy)
        if prefix in grp_readVar:
            if prefix.startswith("MET"):
                grpcls_alt, brMapMap = _makeAltClassAndMaps(
                        grpNm, grp_dict, getJetMETVarName_postproc,
                        systName="jet", nomName="nom", exclVars=("raw",),
                        attCls=altProxy, altBases=(AltLeafGroupProxy,)
                        )
            varsProxy  = AltLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt)
            setTreeAtt(f"_{grpNm}", varsProxy)
            setTreeAtt(grpNm, varsProxy["nomWithSyst"])
            logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in varsProxy["nomWithSyst"].brMap.values())))))

        for vari in systVariations:
            if vari.appliesTo(grpNm):
                grpcls_alt, brMapMap = _makeAltClassAndMaps(
                        grpNm, grp_dict, vari.getVarName,
                        systName=vari.systName, nomName=(vari.origName if vari.isCalc and vari.origName else vari.nomName), exclVars=vari.exclVars,
                        attCls=altProxy, altBases=(AltLeafGroupProxy,)
                        )
                withSyst = "nomWithSyst"
                if vari.isCalc:
                    varsProxy = CalcLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt, withSystName=withSyst)
                else:
                    varsProxy  = AltLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt)
                    logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in varsProxy[withSyst].brMap.values())))))
                setTreeAtt(f"_{grpNm}", varsProxy)
                setTreeAtt(grpNm, varsProxy[withSyst])


    ## SOA, nanoAOD style (LeafCount, shared)
    containerGroupCounts = ("nElectron", "nFatJet", "nIsoTrack", "nJet", "nMuon", "nOtherPV", "nPhoton", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj", "nCorrT1METJet")
    containerGroupCounts_Gen = ("nGenDressedLepton", "nGenJet", "nGenJetAK8", "nGenPart", "nGenVisTau", "nSubGenJetAK8")
    ## check which are there, and for which we need to read variations
    cnt_found = []
    cnt_readVar = []
    for sizeNm in (chain(containerGroupCounts, containerGroupCounts_Gen) if isMC else containerGroupCounts):
        if sizeNm not in allTreeLeafs:
            logger.warning("{0} is not a branch in the tree - skipping collection".format(sizeNm))
        else:
            cnt_found.append(sizeNm)
            if sizeNm == "nJet" and "Jet_pt_nom" in allTreeLeafs:
                cnt_readVar.append(sizeNm)
            if sizeNm == "nMuon" and "Muon_corrected_pt" in allTreeLeafs:
                cnt_readVar.append(sizeNm)

    for sizeNm in cnt_found:
        grpNm = sizeNm[1:]
        prefix = "{0}_".format(grpNm)
        itm_dict = {
            "__doc__" : "{0} proxy class".format(grpNm)
            }
        itm_lvs = set(lvNm for lvNm,lv in allTreeLeafs.items() if lvNm.startswith(prefix) and lv.GetLeafCount().GetName() == sizeNm)
        sizeOp = GetColumn(allTreeLeafs[sizeNm].GetTypeName(), sizeNm)
        for lvNm in itm_lvs:
            lvNm_short = lvNm[len(prefix):]
            col = GetArrayColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm, sizeOp).result
            if "Idx" not in lvNm:
                itm_dict[lvNm_short] = itemProxy(col)
            else:
                coll,i = lvNm_short.split("Idx")
                collPrefix = coll[0].capitalize()+coll[1:]
                collGetter = (
                        GetItemRefCollection_toVar("_{0}".format(collPrefix)) if collPrefix in cnt_readVar or any(vari.appliesTo(collPrefix) for vari in systVariations) else
                        GetItemRefCollection(collPrefix))
                tree_children.append(collGetter)
                itm_dict["".join((coll,i))] = itemRefProxy(col, collGetter)
        ## create p4 branches (naive, but will be reused for variation case)
        if f"{prefix}pt" in itm_lvs and f"{prefix}phi" in itm_lvs:
            itm_dict["p4"] = addP4ToObj(prefix, itm_lvs)
        itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
        ## default collection proxy, replaced below if needed
        coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
        setTreeAtt(grpNm, coll_orig)
        ## insert variations using kinematic calculator, from branches, or not
        if sizeNm in cnt_readVar:
            if sizeNm == "nJet":
                altItemType, brMapMap = _makeAltClassAndMaps(
                        grpNm, itm_dict, getJetMETVarName_postproc,
                        systName="jet", nomName="nom", exclVars=("raw",),
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,)
                        )

            ## add _Jet which holds the variations (not syst-aware), and Jet which is the nominal, with systematics variations (defined just bove)
            varsProxy = AltCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType)
            logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in varsProxy["nomWithSyst"].brMap.values())))))
            setTreeAtt(f"_{grpNm}", varsProxy)
            setTreeAtt(grpNm, varsProxy["nomWithSyst"])

        ## new style generic solution - TODO: implement one case (for each of calc and read)
        for vari in systVariations:
            if vari.appliesTo(grpNm):
                altItemType, brMapMap = _makeAltClassAndMaps(
                        grpNm, itm_dict, vari.getVarName,
                        systName=vari.systName, nomName=(vari.origName if vari.isCalc and vari.origName else vari.nomName), exclVars=vari.exclVars,
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,)
                        )
                withSyst = "nomWithSyst"
                if vari.isCalc:
                    varsProxy = CalcCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType, withSystName=withSyst)
                else:
                    varsProxy = AltCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType)
                    logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in varsProxy[withSyst].brMap.values())))))
                setTreeAtt(f"_{grpNm}", varsProxy)
                setTreeAtt(grpNm, varsProxy[withSyst])

        for lvNm in itm_lvs:
            del tree_dict[lvNm]
        del tree_dict[sizeNm] ## go through op.rng_len

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)
    for pc in tree_children:
        pc._parent = treeProxy

    return treeProxy
