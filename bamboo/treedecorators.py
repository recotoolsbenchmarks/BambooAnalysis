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

class SetAsParentOfAtt(object):
    ## callback to set the tree as _parent of the groups
    def __init__(self, name):
        self.name = name
    def __call__(self, inst):
        getattr(inst, self.name)._parent = inst
class SetAsParent(object):
    def __init__(self, obj):
        self.obj = obj
    def __call__(self, inst):
        self.obj._parent = inst

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
def getJetMETVarName_postproc(nm):
    if nm.split("_")[0] in ("pt", "eta", "phi", "mass") and len(nm.split("_")) == 2:
        return tuple(nm.split("_"))
    elif nm.split("_")[0] == "btagSF":
        if "shape" not in nm.split("_"):
            if "_" not in nm:
                return "btagSF", "nom"
            else: ## btag SF up and down
                return "btagSF", "{0}{1}".format(*nm.split("_"))
        else:
            if nm == "btagSF_shape":
                return "btagSF_shape", "nom"
            elif nm.startswith("btagSF_shape_") and nm.endswith("_jes"):
                return "btagSF_shape", "jesTotal{0}".format(nm.split("_")[2])
            else:
                return "btagSF_shape", "{0}_{1}_{3}{2}".format(*nm.split("_"))
def getPUWeightVarName_postproc(nm):
    if nm.startswith("puWeight"):
        if nm == "puWeight":
            return "puWeight", "nom"
        else:
            return "puWeight", nm

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
        test = getVarName(nm)
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
        SystAltColumnOp(
            getCol(vAtts[nomName]), systName,
            dict((var, getCol(vop)) for var,vop in vAtts.items() if var not in exclVars),
            valid=[ var for var in allVars if var in vAtts ],
            ))
        for attNm,vAtts in var_atts.items())
    return cls_alt, brMapMap

def decorateNanoAOD(aTree, description=None, isMC=False, addCalculators=None):
    """ Decorate a CMS NanoAOD Events tree

    Variation branches following the NanoAODTools conventions (e.g. Jet_pt_nom)
    are automatically used (but calculators for the same collection take
    precendence, if requested).

    :param aTree: TTree to decorate
    :param description: reserved for future version/flavour-dependence
    :param isMC: simulation or not
    :param addCalculators: list of length-branch names for whose collections proxies for a variations calculators should be added
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

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    ## NOTE first attempt: fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm) if not lv.GetLeafCount()
        else GetArrayColumn(lv.GetTypeName(), lvNm, GetColumn(allTreeLeafs[lv.GetLeafCount().GetName()].GetTypeName(), lv.GetLeafCount().GetName()))
        )) for lvNm,lv in allTreeLeafs.items()))
    tree_postconstr = []
    addSetParentToPostConstr = partial(lambda act,obj : act.append(SetAsParent(obj)), tree_postconstr)

    ## weights with variations
    if "puWeightUp" in allTreeLeafs:
        brMap = {}
        toRem = []
        for nm,nmAtt in tree_dict.items():
            test = getPUWeightVarName_postproc(nm)
            if test is not None:
                attNm,varNm = test
                brMap[normVarName(varNm)] = nmAtt.op
                toRem.append(nm)
        for nm in toRem:
            del tree_dict[nm]
        brMap["nomWithSyst"] = SystAltColumnOp(brMap["nom"], "puWeight",
            dict(brMap), valid=[ var for var in brMap.keys() if var != "nom"]
            )
        varsProxy = AltLeafVariations(None, brMap, typeName=brMap["nom"].typeName)
        addSetParentToPostConstr(varsProxy)
        tree_dict["_puWeight"] = varsProxy
        nomSystProxy = varsProxy["nomWithSyst"]
        tree_dict["puWeight"] = nomSystProxy
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
            if prefix.startswith("MET") and "{0}_pt_nom".format(prefix) in allTreeLeafs:
                grp_readVar.append(prefix)
    for prefix in grp_found:
        grpNm = prefix.rstrip("_")
        grp_dict = {
            "__doc__" : "{0} leaf group proxy class".format(grpNm)
            }
        grp_lvNms = set(lvNm for lvNm in allTreeLeafs.keys() if lvNm.startswith(prefix))
        grp_dict.update(dict((lvNm[len(prefix):], proxy(GetColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm))) for lvNm in grp_lvNms))
        grpcls = type("{0}LeafGroupProxy".format(grpNm), (LeafGroupProxy,), grp_dict)
        for lvNm in grp_lvNms:
            del tree_dict[lvNm]
        grp_proxy = grpcls(grpNm, None)
        addSetParentToPostConstr(grp_proxy)
        if prefix.startswith("MET") and addCalculators and prefix.rstrip("_") in addCalculators:
            grpcls_alt, brMapMap = _makeAltClassAndMaps(
                    grpNm, grp_dict, getMETVarName_calc,
                    systName="jet", nomName="raw", exclVars=("raw",),
                    attCls=altProxy, altBases=(AltLeafGroupProxy,)
                    )
            grpNm = "_{0}".format(grpNm)
            varsProxy  = CalcLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt, withSystName="nomWithSyst")
            addSetParentToPostConstr(varsProxy)
            tree_dict[grpNm] = varsProxy
            nomSystProxy = varsProxy["nomWithSyst"]
            addSetParentToPostConstr(nomSystProxy)
            tree_dict[grpNm[1:]] = nomSystProxy
        elif prefix in grp_readVar:
            if prefix.startswith("MET"):
                grpcls_alt, brMapMap = _makeAltClassAndMaps(
                        grpNm, grp_dict, getJetMETVarName_postproc,
                        systName="jet", nomName="nom", exclVars=("raw",),
                        attCls=altProxy, altBases=(AltLeafGroupProxy,)
                        )
            grpNm = "_{0}".format(grpNm)
            varsProxy  = AltLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt)
            addSetParentToPostConstr(varsProxy)
            tree_dict[grpNm] = varsProxy
            nomSystProxy = varsProxy["nomWithSyst"]
            addSetParentToPostConstr(nomSystProxy)
            tree_dict[grpNm[1:]] = nomSystProxy
            logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in nomSystProxy.brMap.values())))))
        else:
            tree_dict[grpNm] = grp_proxy

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
            if not ( addCalculators and sizeNm in addCalculators ):
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
                        GetItemRefCollection_toVar("_{0}".format(collPrefix)) if collPrefix in cnt_readVar or ( addCalculators and collPrefix in addCalculators ) else
                        GetItemRefCollection(collPrefix))
                addSetParentToPostConstr(collGetter)
                itm_dict["".join((coll,i))] = itemRefProxy(col, collGetter)
        ## create p4 branches (naive, but will be reused for variation case)
        p4AttNames = ("pt", "eta", "phi", "mass")
        if all(("".join((prefix, att)) in itm_lvs) for att in p4AttNames):
            itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst.pt, inst.eta, inst.phi, inst.mass)).result)
        itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
        ## insert variations using kinematic calculator, from branches, or not
        if addCalculators and sizeNm in addCalculators:
            coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            addSetParentToPostConstr(coll_orig)
            if sizeNm == "nJet":
                altItemType, brMapMap = _makeAltClassAndMaps(
                        grpNm, itm_dict, getJetEnergyVarName_calc,
                        systName="jet", nomName="raw", exclVars=("raw",),
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,)
                        )
            elif sizeNm == "nMuon":
                altItemType, brMapMap = _makeAltClassAndMaps(
                        grpNm, itm_dict, getMuonMomentumVarName_calc,
                        systName="muon", nomName="raw", exclVars=("raw",),
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,)
                        )
            else:
                raise RuntimeError("Adding a calculator for collection {0} is not supported (yet)".format(sizeNm))
            grpNm = "_{0}".format(grpNm) ## add variations as '_Muon'/'_Jet', nominal as 'Muon', 'Jet'
            varsProxy = CalcCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType, withSystName="nomWithSyst")
            tree_dict[grpNm] = varsProxy
            nomSystProxy = varsProxy["nomWithSyst"]
            addSetParentToPostConstr(nomSystProxy)
            tree_dict[grpNm[1:]] = nomSystProxy
        elif sizeNm in cnt_readVar:
            coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            addSetParentToPostConstr(coll_orig)
            if sizeNm == "nJet":
                altItemType, brMapMap = _makeAltClassAndMaps(
                        grpNm, itm_dict, getJetMETVarName_postproc,
                        systName="jet", nomName="nom", exclVars=("raw",),
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,)
                        )

            ## add _Jet which holds the variations (not syst-aware), and Jet which is the nominal, with systematics variations (defined just bove)
            grpNm = "_{0}".format(grpNm)
            tree_dict[grpNm] = AltCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType)
            nomSystProxy = tree_dict[grpNm]["nomWithSyst"]
            addSetParentToPostConstr(nomSystProxy)
            tree_dict[grpNm[1:]] = nomSystProxy
            logger.debug("{0} variations read from branches: {1}".format(grpNm, list(set(chain.from_iterable(op.variations for op in nomSystProxy.brMap.values())))))
        else:
            tree_dict[grpNm] = ContainerGroupProxy(prefix, None, sizeOp, itmcls)

        for lvNm in itm_lvs:
            del tree_dict[lvNm]
        del tree_dict[sizeNm] ## go through op.rng_len
        tree_postconstr.append(SetAsParentOfAtt(grpNm))

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)

    for pc in tree_postconstr:
        pc(treeProxy)

    return treeProxy
