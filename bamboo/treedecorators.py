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

def decorateNanoAOD(aTree, description=None, isMC=False):
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
    for sizeNm in (chain(containerGroupCounts, containerGroupCounts_Gen) if isMC else containerGroupCounts):
        if sizeNm not in allTreeLeafs:
            logger.warning("{} is not a branch in the tree!".format(sizeNm))
            continue
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
                        GetItemRefCollection_toVar("_{0}".format(collPrefix)) if collPrefix in ("Jet", "Muon") else
                        GetItemRefCollection(collPrefix))
                addSetParentToPostConstr(collGetter)
                itm_dict["".join((coll,i))] = itemRefProxy(col, collGetter)
        p4AttNames = ("pt", "eta", "phi", "mass")
        if all(("".join((prefix, att)) in itm_lvs) for att in p4AttNames):
            if sizeNm not in ("nJet", "nMuon"):
                itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst.pt, inst.eta, inst.phi, inst.mass)).result)
            else:
                for att in p4AttNames:
                    itm_dict["_{0}".format(att)] = itm_dict[att]
                itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst._pt, inst._eta, inst._phi, inst._mass)).result) # note not too efficient, but only for "original" jet collection ("nominal" should be the default)
        itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
        if sizeNm in ("nJet", "nMuon"):
            coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            addSetParentToPostConstr(coll_orig)
            itm_dict_var = _translate_to_var(itm_dict, toSkip=("p4",), addToPostConstr=addSetParentToPostConstr)
            itm_dict_var["p4"] = funProxy(lambda inst : inst._parent._parent.result.momenta()[inst._idx])
            varItemType = type("Var{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict_var)
            # construct arguments for correction/variation calculator
            anObj = coll_orig[0]
            args = [ comp.op.arg for comp in (anObj.pt, anObj.eta, anObj.phi, anObj.mass) ]
            nameMap = None
            import copy
            if sizeNm == "nJet":
                if isMC:
                    ## not the most elegant solution, but we can't make assumptions on the relative order with respect to and MET
                    genkinj = copy.deepcopy(args)
                    for it in genkinj:
                        it.name = it.name.replace("Jet", "GenJet")
                else:
                    genkinj = list(repeat(ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}"), 4))
                args += [ comp.op.arg for comp in (anObj.rawFactor, anObj.area) ]
                args += [ GetColumn("Float_t", nm) for nm in ("fixedGridRhoFastjetAll", "MET_phi", "MET_pt", "MET_sumEt") ]
                args += genkinj
            elif sizeNm == "nMuon":
                args += [ comp.op.arg for comp in (anObj.charge, anObj.nTrackerLayers) ]
                if isMC:
                    genpi = copy.deepcopy(anObj.nTrackerLayers.op.arg)
                    genpi.name = genpi.name.replace("nTrackerLayers", "genPartIdx")
                    args.append(genpi)
                    genpt = copy.deepcopy(anObj.pt.op.arg)
                    genpt.name = genpt.name.replace("Muon", "GenPart")
                    args.append(genpt)
                else:
                    args.append(ExtVar("ROOT::VecOps::RVec<Int_t>", "ROOT::VecOps::RVec<Int_t>{}"))
                    args.append(ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}"))
            nameMap={"nominal": grpNm}
            grpNm = "_{0}".format(grpNm) ## add variations as '_Muon'/'_Jet', nominal as 'Muon', 'Jet'
            tree_dict[grpNm] = Variations(None, coll_orig, args, varItemType=varItemType, nameMap=nameMap)
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
