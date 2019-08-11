"""
Main user-visible module of tree operations, separate from the operation definitions
and proxy classes
User-facing module
"""

from functools import partial

from .treeproxies import *
from . import treefunctions as op

logger = logging.getLogger(__name__)

def allLeafs(branch):
    """
    Recursively collect TTree leaves (TLeaf and TBranchElement)
    """
    for br in branch.GetListOfBranches():
        from cppyy import gbl
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

class varItemProxy(object):
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op[inst._parent._parent.result.indices()[inst._idx]]
class varItemRefProxy(object):
    def __init__(self, op, getTarget):
        self.op = op
        self.getTarget = getTarget
    def __get__(self, inst, cls):
        return self.getTarget(inst)[self.op[inst._parent._parent.result.indices()[inst._idx]]]

class altItemProxy(object):
    def __init__(self, name, op):
        self.name = name
        self.op = op
    def __get__(self, inst, cls):
        return inst._parent.brMap.get(self.name, op)[inst._idx]

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
            from cppyy import gbl
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

    def _translate_to_var(itm_dict, toSkip=None, addToPostConstr=None):
        ## translate attributes NOTE there may be a more elegant solution to this,
        ## like moving the logic into ContainerGroupItemProxy
        ## (similar to ListBase: know about the base collection and the final index there as well)
        itm_dict_var = dict()
        for nm,itmAtt in itm_dict.items():
            if isinstance(itmAtt, itemProxy): ## regular branch, change just index
                itm_dict_var[nm] = varItemProxy(itmAtt.op)
            elif isinstance(itmAtt, itemRefProxy):
                itm_dict_var[nm] = varItemRefProxy(itmAtt.op, itmAtt.getTarget)
            elif isinstance(itmAtt, funProxy) and toSkip is not None and nm in toSkip:
                pass ## ok to skip
            elif isinstance(itmAtt, str): ## __doc__
                itm_dict_var[nm] = itmAtt
            else:
                raise RuntimeError("Cannot translate attribute {0} for variations yet".format(nm))
        return itm_dict_var

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    ## NOTE first attempt: fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm))) for lvNm,lv in allTreeLeafs.items()))
    tree_postconstr = []
    addSetParentToPostConstr = partial(lambda act,obj : act.append(SetAsParent(obj)), tree_postconstr)
    simpleGroupPrefixes = ("CaloMET_", "ChsMET_", "MET_", "PV_", "PuppiMET_", "RawMET_", "TkMET_", "Flag_", "HLT_") ## TODO get this from description?
    for prefix in (chain(simpleGroupPrefixes, ("GenMET_", "Generator_", "LHE_",)) if isMC else simpleGroupPrefixes):
        grpNm = prefix.rstrip("_")
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
    ## SOA, nanoAOD style (LeafCount, shared)
    containerGroupCounts = ("nElectron", "nFatJet", "nIsoTrack", "nJet", "nMuon", "nOtherPV", "nPhoton", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj")
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
                    logger.debug("Will read Jet variations from branches")
                    cnt_readVar.append(sizeNm)
                if sizeNm == "nMuon" and "Muon_corrected_pt" in allTreeLeafs:
                    logger.debug("Will read Muon variations from branches")
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
            itm_dict_var = _translate_to_var(itm_dict, toSkip=("p4",), addToPostConstr=addSetParentToPostConstr)
            itm_dict_var["p4"] = funProxy(lambda inst : inst._parent._parent.result.momenta()[inst._idx])
            varItemType = type("Var{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict_var)
            # construct arguments for correction/variation calculator
            nameMap={"nominal": grpNm}
            grpNm = "_{0}".format(grpNm) ## add variations as '_Muon'/'_Jet', nominal as 'Muon', 'Jet'
            tree_dict[grpNm] = CalcVariations(None, coll_orig, varItemType=varItemType, nameMap=nameMap)
        elif sizeNm in cnt_readVar:
            coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            addSetParentToPostConstr(coll_orig)
            if sizeNm == "nJet":
                ## collect ops of kinematic variables that change (nominal as well as varied)
                pt_atts = dict((nm, att.op) for nm,att in itm_dict.items() if nm.startswith("pt_"))
                mass_atts  = dict((nm, att.op) for nm,att in itm_dict.items() if nm.startswith("mass_"))
                ## redirect in altItemproxy
                itm_dict_alt = dict(itm_dict)
                itm_dict_alt["pt"] = altItemProxy("pt", itm_dict["pt"].op)
                itm_dict_alt["mass"] = altItemProxy("mass", itm_dict["mass"].op)
                for nm in chain(pt_atts.keys(), mass_atts.keys()):
                    del itm_dict_alt[nm]
                pt_atts = dict((normVarName(nm[3:]), att) for nm,att in pt_atts.items())
                mass_atts = dict((normVarName(nm[5:]), att) for nm,att in mass_atts.items())

                altItemType = type("Alt{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict_alt)
                ## construct the map of maps of redirections { variation : { varName : op } }
                brMapMap = {}
                for var,vop in pt_atts.items():
                    if var not in brMapMap:
                        brMapMap[var] = {}
                    brMapMap[var]["pt"] = vop
                for var,vop in mass_atts.items():
                    if var not in brMapMap:
                        brMapMap[var] = {}
                    brMapMap[var]["mass"] = vop
                ## nominal: with systematic variations (all are valid, but not all need to modify)
                allVars = list(k for k in brMapMap.keys() if k not in ("raw", "nom"))
                brMapMap["nomWithSyst"] = {
                    "pt" : SystAltColumnOp(pt_atts["nom"].op, "jet", dict((var, vop.op.name) for var,vop in pt_atts.items() if var != "raw"), valid=allVars).result,
                    "mass"  : SystAltColumnOp(mass_atts["nom"].op, "jet", dict((var, vop.op.name) for var,vop in mass_atts.items() if var != "raw"), valid=allVars).result
                    }
                ## add _Jet which holds the variations (not syst-aware), and Jet which is the nominal, with systematics variations (defined just bove)
                grpNm = "_{0}".format(grpNm)
                tree_dict[grpNm] = AltVariations(None, coll_orig, brMapMap, altItemType=altItemType)
                nomSystProxy = tree_dict[grpNm]["nomWithSyst"]
                tree_postconstr.append(SetAsParent(nomSystProxy))
                tree_dict[grpNm[1:]] = nomSystProxy
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
