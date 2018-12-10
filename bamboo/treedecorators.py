"""
Main user-visible module of tree operations, separate from the operation definitions
and proxy classes
User-facing module
"""

from .treeproxies import *
from . import treefunctions as op

def allLeafs(branch, split=False): ## study splitlevel in tree as well
    """
    Recursively collect TTree leaves (TLeaf and TBranchElement)
    """
    if split:
        raise NotImplementedError("Going into TBranchElement not implemented yet")
    for lv in branch.GetListOfLeaves():
        yield lv
    for br in branch.GetListOfBranches():
        import ROOT
        if isinstance(br, ROOT.TBranchElement):
            yield br
        else:
            for lv in allLeafs(br):
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
            import ROOT
            for lvNm, lv in allTreeLeafs.items():
                if lvNm.startswith(prefix):
                    if isinstance(lv, ROOT.TBranchElement):
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
                        print("DEBUG: detected a split LorentzVector array {}".format(cntLvNm.rstrip("_")))
                        atts = [ aLvNm[len(cntLvNm)+len("fCoordinates.f"):] for aLvNm in arrLvs ]
                        print(" ".join(aLvNm[len(cntLvNm)+len("fCoordinates.f"):] for aLvNm in arrLvs))
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
                        print("Warning: Detected a split object that's not a LorentzVector. Attributes are: {}".format(", ".join(aLvNm[len(cntLvNm):] for aLvNm in arrLvs)))
                else:
                    for arrLv in arrLvs:
                        lvNm_short = arrLv[len(prefix):]
                        col = GetArrayColumn(allTreeLeafs[arrLv].GetTypeName(), arrLv, sizeOp).result
                        itm_dict[lvNm_short] = itemProxy(col)
                        del tree_dict[arrLv]
                    del tree_dict[cntLvNm]
            aLvNm = next(lv for lv in itm_lvs_vec if "{0}_".format(lv) not in itm_lvs_vec)
            sizeOp = adaptArg(op.rng_len(GetColumn(allTreeLeafs[aLvNm].GetTypeName(), aLvNm).result)) ## abuse a bit
            for lvNm in itm_lvs_vec:
                lvNm_short = lvNm[len(prefix):]
                col = GetColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm).result
                itm_dict[lvNm_short] = itemProxy(col)
                del tree_dict[lvNm]
            itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
            tree_dict[grpNm] = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            tree_postconstr.append(SetAsParentOfAtt(grpNm))
            print("Processed group {0}, remaining description: {1}".format(grpNm, desc_rem))
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

    print(tree_dict)

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)

    for pc in tree_postconstr:
        pc(treeProxy)

    return treeProxy

def decorateNanoAOD(aTree, description=None):
    if description is None:
        description = dict()

    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    ## NOTE first attempt: fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm))) for lvNm,lv in allTreeLeafs.items()))
    tree_postconstr = []
    simpleGroupPrefixes = ("CaloMET_", "Flag_", "HLT_", "RawMET_", "PuppiMET_", "MET_", "TkMET_", "PV_") ## TODO get this from description
    for prefix in simpleGroupPrefixes:
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
    containerGroupCounts = ("nElectron", "nMuon", "nTau", "nPhoton", "nJet", "nFatJet", "nSubJet", "nTrigObj", "nSV", "nOtherPV", "nSoftActivityJet")
    for sizeNm in containerGroupCounts:
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
                itm_dict["".join((coll,i))] = itemRefProxy(col, (lambda colName : (lambda me : getattr(me._parent._parent, colName))) (coll.capitalize()))
        p4AttNames = ("pt", "eta", "phi", "mass")
        if all(("".join((prefix, att)) in itm_lvs) for att in p4AttNames):
            if sizeNm != "nJet":
                itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst.pt, inst.eta, inst.phi, inst.mass)).result)
            else:
                for att in p4AttNames:
                    itm_dict["_{0}".format(att)] = itm_dict[att]
                itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst._pt, inst._eta, inst._phi, inst._mass)).result) ## TODO replace by get from the right vector
        itmcls = type("{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict)
        if sizeNm != "nJet":
            tree_dict[grpNm] = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
        else:
            jets_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
            tree_postconstr.append(SetAsParent(jets_orig))
            ## translate attributes NOTE there may be a more elegant solution to this,
            ## like moving the logic into ContainerGroupItemProxy
            ## (similar to ListBase: know about the base collection and the final index there as well)
            itm_dict_var = dict()
            for nm,itmAtt in itm_dict.items():
                if isinstance(itmAtt, itemProxy): ## regular branch, change just index
                    itm_dict_var[nm] = varItemProxy(itmAtt.op)
                elif isinstance(itmAtt, itemRefProxy):
                    itm_dict_var[nm] = varItemRefProxy(itmAtt.op, itmAtt.getTarget)
                elif isinstance(itmAtt, funProxy) and ( nm == "p4" ):
                    pass ## ok to skip
                elif isinstance(itmAtt, str): ## __doc__
                    itm_dict_var[nm] = itmAtt
                else:
                    raise RuntimeError("Cannot translate attribute {0} for jet variations yet".format(nm))
            itm_dict_var["p4"] = funProxy(lambda inst : inst._parent._parent.momenta()[inst._idx]) ## container parent is the modifier op
            varItemType = type("Var{0}GroupItemProxy".format(grpNm), (ContainerGroupItemProxy,), itm_dict_var)

            tree_dict[grpNm] = Variations(None, jets_orig, varItemType=varItemType)

        for lvNm in itm_lvs:
            del tree_dict[lvNm]
        del tree_dict[sizeNm] ## go through op.rng_len
        tree_postconstr.append(SetAsParentOfAtt(grpNm))

    TreeProxy = type("{0}Proxy".format(aTree.GetName()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)

    for pc in tree_postconstr:
        pc(treeProxy)

    return treeProxy