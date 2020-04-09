"""
Expressions are constructed by executing python code on decorated versions of
decorated trees. The :py:mod:`bamboo.treedecorators` module contains helper
methods to do so for commonly used formats, e.g. :py:func:`~.decorateNanoAOD`
for CMS NanoAOD.
"""

from collections import defaultdict
from functools import partial

from .treeproxies import *
from . import treefunctions as op

import logging
logger = logging.getLogger(__name__)

def allLeafs(branch):
    # Recursively collect TTree leaves (TLeaf and TBranchElement)
    for br in branch.GetListOfBranches():
        from .root import gbl
        if isinstance(br, gbl.TBranchElement):
            yield br
        else:
            for lv in br.GetListOfLeaves():
                yield lv

def normVarName(varName):
    # Normalize variation name: if ending in up or down, make sure this part has no capitals (for plotIt)
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

def _makeAltClassAndMaps(name, dict_orig, vari, getCol=lambda op : op, attCls=None, altBases=None): ## internal, leaf/group/collection
    ## vari.getVarName should return the variable and variation name (nomName for the nominal one)
    ## if this is a systematic variation branch - otherwise None
    dict_alt = dict(dict_orig)
    ## collect ops of kinematic variables that change (nominal as well as varied)
    var_atts = defaultdict(dict)
    for nm, nmAtt in dict_orig.items():
        test = vari.getVarName(nm, collgrpname=name)
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
    nomName = vari.origName if vari.isCalc else vari.nomName(name)
    exclVars = list(normVarName(var) for var in vari.exclVars(name))
    allVars = list(k for k in brMapMap.keys() if k not in exclVars and k != nomName)
    brMapMap["nomWithSyst"] = dict((attNm,
        SystAltOp(
            getCol(vAtts[nomName]), vari.systName,
            dict((var, getCol(vop)) for var,vop in vAtts.items() if var not in exclVars and var != nomName),
            valid=tuple(var for var in allVars if var in vAtts),
            ))
        for attNm,vAtts in var_atts.items())
    return cls_alt, brMapMap

## Helper classes
class NanoSystematicVarSpec:
    """ Interface for classes that specify how to incorporate systematics or on-the-fly corrections in the decorated tree

    See :py:class:`~.NanoAODDescription` and :py:func:`~.decorateNanoAOD`
    """
    def __init__(self, systName, nomName=None, origName=None, exclVars=None, isCalc=False):
        """ Base class constructor

        :param systName: name of the systematic source, for automatic systematics (if applicable)
        :param nomName: nominal variation name (for non-calculated variations; can be customised in subclasses by overriding the :py:meth:`~.nomName` method)
        :param origName: original variation name (for calculated variations)
        :param exclVars: variations that are found but should not be used in automatic systematics (can be customised through :py:meth:`~.nomName`)
        :param isCalc: boolean indicating whether variations are calculated or read from alternative branches
        """
        self.systName = systName
        self._nomName = nomName
        self.origName = origName
        self._exclVars = exclVars if exclVars is not None else tuple()
        self.isCalc = isCalc
    def appliesTo(self, name):
        """ Return true if this systematic variation requires action for this variable, group, or collection """
        return False
    def getVarName(self, branchName, collgrpname=None):
        """ Get the variable name and variation corresponding to an (unprefixed, in case of groups or collections) branch name """
        pass
    def nomName(self, name):
        """ Nominal systematic variation name for a group/collection """
        return self._nomName
    def exclVars(self, name):
        """ Systematic variations to exclude for a group/collection """
        return self._exclVars

class ReadVariableVarWithSuffix(NanoSystematicVarSpec):
    """ Read variations of a single branch from branches with the same name with a suffix """
    def __init__(self, commonName, sep="_", systName=None, nomName=None, exclVars=None):
        super(ReadVariableVarWithSuffix, self).__init__(
                (systName if systName is not None else commonName),
                nomName=nomName, exclVars=exclVars, isCalc=False)
        self.prefix = commonName
        self.sep = sep
    def appliesTo(self, name):
        """ True if name starts with the prefix """
        return name.startswith(self.prefix)
    def getVarName(self, branchName, collgrpname=None):
        """ Split into prefix and variation (if present, else nominal) """
        variNm = normVarName(branchName[len(self.prefix):].lstrip(self.sep))
        return self.prefix, variNm if variNm else self.nomName

nanoPUWeightVar = ReadVariableVarWithSuffix("puWeight")

class ReadJetMETVar(NanoSystematicVarSpec):
    """
    Read jet and MET kinematic variations from different branches for automatic systematic uncertainties

    :param jetsName: jet collection prefix (e.g. ``"Jet"``)
    :param metName: MET prefix (e.g. ``"MET"``)
    :param systName: systematic group name (backend optimisation detail, ``"jet"`` now)
    :param jetsNomName: name of the nominal jet variation (``"nom"`` by default)
    :param jetsOrigName: name of the original jet variation (``"raw"`` by default)
    :param metNomName: name of the nominal jet variation (``"nom"`` by default)
    :param metOrigName: name of the original jet variation (``"raw"`` by default)
    :param jetsExclVars: jet variations that are present but should be ignored (if not specified, only ``jetsOrigName`` is taken, so if specified this should usually be added explicitly)
    :param metExclVars: MET variations that are present but should be ignored (if not specified, only ``metOrigName`` is taken, so if specified this should usually be added explicitly)
    :param bTaggers: list of b-tagging algorithms, for scale factors stored in a branch
    :param bTagWPs: list of b-tagging working points, for scale factors stored in a branch (``shape`` should be included here, if wanted)

    .. note:: The implementation of automatic systematic variations treats
       "xyzup" and "xyzdown" independently (since this is the most flexible).
       If a source of systematic uncertainty should be excluded, both the "up"
       and "down" variation should then be added to the list of variations to
       exclude (``jetsExclVars`` or ``metExclVars``).
    """
    def __init__(self, jetsName, metName, systName="jet", jetsNomName="nom", jetsOrigName="raw", metNomName="jer", metOrigName="raw", jetsExclVars=None, metExclVars=None, bTaggers=None, bTagWPs=None):
        super(ReadJetMETVar, self).__init__(
                systName, nomName=jetsNomName, origName=jetsOrigName, isCalc=False,
                exclVars=(jetsExclVars if jetsExclVars is not None else (jetsOrigName,)))
        self.jetsName = jetsName
        self.metName = metName
        self.metOrigname = metOrigName
        self.bTaggers = bTaggers if bTaggers is not None else []
        self.bTagWPs = bTagWPs if bTagWPs is not None else []
        self.metNomName = metNomName
        self.metExclVars = (metExclVars if metExclVars is not None else (metOrigName,))
    def appliesTo(self, name):
        return name in (self.jetsName, self.metName)
    def nomName(self, name):
        if name == self.metName:
            return self.metNomName
        else:
            return self._nomName
    def exclVars(self, name):
        if name == self.metName:
            return self.metExclVars
        else:
            return self._exclVars
    def getVarName(self, nm, collgrpname=None):
        if nm.split("_")[0] in ("pt", "eta", "phi", "mass") and len(nm.split("_")) >= 2:
            return (nm.split("_")[0], "_".join(nm.split("_")[1:]))
        elif nm.startswith("btagSF"):
            for tagger in self.bTaggers:
                for wp in self.bTagWPs:
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

nanoReadJetMETVar = ReadJetMETVar("Jet", "MET", bTaggers=["csvv2", "deepcsv", "deepjet", "cmva"], bTagWPs=["L", "M", "T", "shape"])
nanoReadJetMETVar_METFixEE2017 = ReadJetMETVar("Jet", "METFixEE2017", bTaggers=["csvv2", "deepcsv", "deepjet", "cmva"], bTagWPs=["L", "M", "T", "shape"])

class CalcCollection(NanoSystematicVarSpec):
    """ Base class for :py:class:`~.NanoSystematicVarSpec` with a calculator """
    def __init__(self, collName, calcAttrs, systName=None, nomName="nominal", origName="raw", exclVars=None):
        super(CalcCollection, self).__init__(
                (systName if systName is not None else collName.lower()),
                exclVars=(exclVars if exclVars is not None else (origName,)),
                nomName=nomName, origName=origName, isCalc=True)
        self.collName = collName
        self.calcAttrs = calcAttrs
    def appliesTo(self, name):
        """ True if name equals the collection name """
        return name == self.collName
    def getVarName(self, nm, collgrpname=None):
        """ Valid if name is a calculated attribute """
        if nm in self.calcAttrs:
            return nm, self.origName

class CalcJetMETVar(NanoSystematicVarSpec):
    """ :py:class:`~.NanoSystematicVarSpec` for on-the-fly Jet and MET correction and systematic variation calculation """
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

class NanoAODDescription:
    """ Description of the expected NanoAOD structure, and configuration for systematics and corrections

    Essentially, a collection of three containers:
       - :py:attr:`~.collections` a list of collections (by the name of the length leaf)
       - :py:attr:`~.groups` a list of non-collection groups (by prefix, e.g. ``HLT_``)
       - :py:attr:`~.systVariations` a list of :py:class:`~.NanoSystematicVarSpec` instances, to configure reading systematics variations from branches, or calculating them on the fly

    The recommended way to obtain a configuration is from the factory method :py:meth:`~.get`
    """
    ProductionVersions = dict()
    def __init__(self, groups=None, collections=None, systVariations=None):
        self.groups = list(groups) if groups is not None else []
        self.collections = list(collections) if collections is not None else []
        self.systVariations = list(systVariations) if systVariations is not None else []
    def clone(self, addGroups=None, removeGroups=None, addCollections=None, removeCollections=None, systVariations=None):
        groups = self.groups
        if removeGroups:
            groups = [ grp for grp in groups if grp not in removeGroups ]
        if addGroups:
            groups = groups + [ grp for grp in addGroups if grp not in groups ]
        collections = self.collections
        if removeCollections:
            collections = [ grp for grp in collections if grp not in removeCollections ]
        if addCollections:
            collections = collections + [ grp for grp in addCollections if grp not in collections ]
        return self.__class__(groups=groups, collections=collections, systVariations=systVariations)
    @staticmethod
    def get(tag, year="2016", isMC=False, addGroups=None, removeGroups=None, addCollections=None, removeCollections=None, systVariations=None):
        """ Create a suitable NanoAODDescription instance based on a production version

        A production version is defined by a tag, data-taking year, and a flag
        to distinguish data from simulation.
        Any number of groups or collections can be added or removed from this.
        The ``systVariations`` option

        :Example:

        >>> decorateNanoAOD(tree, NanoAODDescription.get("v5", year="2016", isMC=True, systVariations=[nanoRochesterCalc, nanoJetMETCalc]))
        >>> decorateNanoAOD(tree, NanoAODDescription.get("v5", year="2017", isMC=True, systVariations=[nanoPUWeightVar, nanoReadJetMETVar_METFixEE2017]))

        :param tag: production version (e.g. "v5")
        :param year: data-taking year
        :param isMC: simulation or not
        :param addGroups: (optional) list of groups of leaves to add (e.g. ``["L1_", "HLT_"]``, if not present)
        :param removeGroups: (optional) list of groups of leaves to remove (e.g. ``["L1_"]``, if skimmed)
        :param addCollections: (optional) list of containers to add (e.g. ``["nMyJets"]``)
        :param removeCollections: (optional) list of containers to remove (e.g. ``["nPhoton", "nTau"]``)
        :param systVariations: list of correction or systematic variation on-the-fly calculators or configurations to add (:py:class:`~.NanoSystematicVarSpec` instances)

        See also :py:func:`~.decorateNanoAOD`
        """
        return NanoAODDescription.ProductionVersions[(tag, year, isMC)].clone(addGroups=addGroups, removeGroups=removeGroups, addCollections=addCollections, removeCollections=removeCollections, systVariations=systVariations)

_ndpv = NanoAODDescription.ProductionVersions
_ndpv[("v5", "2016", False)] = NanoAODDescription(
        groups=["CaloMET_", "ChsMET_", "MET_", "PV_", "PuppiMET_", "RawMET_", "TkMET_", "Flag_", "HLT_", "L1_"],
        collections=["nElectron", "nFatJet", "nIsoTrack", "nJet", "nMuon", "nOtherPV", "nPhoton", "nSV", "nSoftActivityJet", "nSubJet", "nTau", "nTrigObj", "nCorrT1METJet"])
_ndpv[("v5", "2016", True )] = _ndpv[("v5", "2016", False)].clone(
        addGroups=["GenMET_", "Generator_", "LHE_", "HTXS_"],
        addCollections=["nGenDressedLepton", "nGenJet", "nGenJetAK8", "nGenPart", "nGenVisTau", "nSubGenJetAK8"])
_ndpv[("v5", "2017", False)] = _ndpv[("v5", "2016", False)].clone(addGroups=["METFixEE2017_"])
_ndpv[("v5", "2017", True )] = _ndpv[("v5", "2016", True )].clone(addGroups=["METFixEE2017_"])
_ndpv[("v5", "2018", False)] = _ndpv[("v5", "2016", False)]
_ndpv[("v5", "2018", True )] = _ndpv[("v5", "2016", True )]

def decorateNanoAOD(aTree, description=None):
    """ Decorate a CMS NanoAOD Events tree

    Variation branches following the NanoAODTools conventions (e.g. Jet_pt_nom)
    are automatically used (but calculators for the same collection take
    precendence, if requested).

    :param aTree: TTree to decorate
    :param description: description of the tree format, and configuration for reading or calculating systematic variations and corrections, a :py:class:`~.NanoAODDescription` instance (see also :py:meth:`.NanoAODDescription.get`)
    """
    if description is None:
        raise ValueError('A description is needed to correctly decorate a NanoAOD (but it may be as simple as ``description=NanoAODDescription.get("v5")``)')

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
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm) if not lv.GetLeafCount()
        else GetArrayColumn(lv.GetTypeName(), lvNm, GetColumn(allTreeLeafs[lv.GetLeafCount().GetName()].GetTypeName(), lv.GetLeafCount().GetName()))
        )) for lvNm,lv in allTreeLeafs.items()))
    tree_children = list()
    def setTreeAtt(name, proxy, setParent=True):
        tree_dict[name] = proxy
        if setParent:
            tree_children.append(proxy)

    ## variables with variations
    for vari in description.systVariations:
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
    grp_found = []
    for prefix in description.groups:
        if not any(lvNm.startswith(prefix) for lvNm in allTreeLeafs):
            logger.warning("No branch name starting with {0} in the tree - skipping group".format(prefix))
        else:
            grp_found.append(prefix)
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
        for vari in description.systVariations:
            if vari.appliesTo(grpNm):
                grpcls_alt, brMapMap = _makeAltClassAndMaps(grpNm, grp_dict, vari,
                        attCls=altProxy, altBases=(AltLeafGroupProxy,))
                withSyst = "nomWithSyst"
                if vari.isCalc:
                    varsProxy = CalcLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt, withSystName=withSyst)
                else:
                    varsProxy  = AltLeafGroupVariations(None, grp_proxy, brMapMap, grpcls_alt)
                    allVars = list(set(chain.from_iterable(op.variations for op in varsProxy[withSyst].brMap.values())))
                    if allVars:
                        logger.debug(f"{grpNm} variations read from branches: {allVars}")
                setTreeAtt(f"_{grpNm}", varsProxy)
                setTreeAtt(grpNm, varsProxy[withSyst])

    class NanoAODGenRanges:
        @property
        def parent(self):
            return self.genPartMother
        @property
        def ancestors(self):
            return SelectionProxy(self._parent, Construct("rdfhelpers::gen::ancestors<-1,ROOT::VecOps::RVec<Int_t>>", (op.static_cast("int", self.genPartMother._idx), self.genPartMother._idx.arg)), type(self))

    ## SOA, nanoAOD style (LeafCount, shared)
    cnt_found = []
    for sizeNm in description.collections:
        if sizeNm not in allTreeLeafs:
            logger.warning("{0} is not a branch in the tree - skipping collection".format(sizeNm))
        else:
            cnt_found.append(sizeNm)

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
                        GetItemRefCollection_toVar("_{0}".format(collPrefix)) if any(vari.appliesTo(collPrefix) for vari in description.systVariations) else
                        GetItemRefCollection(collPrefix))
                tree_children.append(collGetter)
                itm_dict["".join((coll,i))] = itemRefProxy(col, collGetter)
        ## create p4 branches (naive, but will be reused for variation case)
        if f"{prefix}pt" in itm_lvs and f"{prefix}phi" in itm_lvs:
            itm_dict["p4"] = addP4ToObj(prefix, itm_lvs)
        itm_bases = [ContainerGroupItemProxy]
        if sizeNm == "nGenPart":
            itm_bases.append(NanoAODGenRanges)
        itmcls = type("{0}GroupItemProxy".format(grpNm), tuple(itm_bases), itm_dict)
        ## default collection proxy, replaced below if needed
        coll_orig = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
        setTreeAtt(grpNm, coll_orig)
        ## insert variations using kinematic calculator, from branches, or not
        for vari in description.systVariations:
            if vari.appliesTo(grpNm):
                altItemType, brMapMap = _makeAltClassAndMaps(grpNm, itm_dict, vari,
                        getCol=(lambda att : att.op), attCls=altItemProxy, altBases=(ContainerGroupItemProxy,))
                withSyst = "nomWithSyst"
                if vari.isCalc:
                    varsProxy = CalcCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType, withSystName=withSyst)
                else:
                    varsProxy = AltCollectionVariations(None, coll_orig, brMapMap, altItemType=altItemType)
                    allVars = list(set(chain.from_iterable(op.variations for op in varsProxy[withSyst].brMap.values())))
                    if allVars:
                        logger.debug(f"{grpNm} variations read from branches: {allVars}")
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

def decorateCMSPhase2SimTree(aTree, isMC=True):
    """ Decorate a flat tree as used for CMS Phase2 physics studies """

    class Phase2GenRanges:
        @property
        def parent(self):
            return self._parent[self.m1]
        @property
        def ancestors(self):
            return SelectionProxy(self._parent, Construct("rdfhelpers::gen::ancestors<-99>", (op.static_cast("int", self.parent._idx), self.parent._idx.arg)), type(self))
        @property
        def children(self):
            return self._parent[self.d1:op.switch(self.d2 > 0, self.d2+1, self.d1)]
        @property
        def descendants(self):
            return SelectionProxy(self._parent, Construct("rdfhelpers::gen::descendants_firstlast<-99,ROOT::VecOps::RVec<Int_t>>", (op.static_cast("int", self._idx), self.d1._parent.arg, self.d2._parent.arg)), type(self))

    known_one = ("event",)
    allTreeLeafs = dict((lv.GetName(), lv) for lv in allLeafs(aTree))
    tree_dict = {"__doc__" : "{0} tree proxy class".format(aTree.GetName())}
    ## fill all, take some out later
    tree_dict.update(dict((lvNm, proxy(GetColumn(lv.GetTypeName(), lvNm))) for lvNm,lv in allTreeLeafs.items()))
    tree_children = []
    cnt_lvs = set(lv.GetLeafCount().GetName() for lv in allTreeLeafs.values() if lv.GetLeafCount())
    for sizeNm in cnt_lvs:
        grpNm = sizeNm.split("_")[0]
        prefix = "{0}_".format(grpNm)
        itm_dict = {
            "__doc__" : "{0} proxy class".format(grpNm)
            }
        itm_lvs = set(lvNm for lvNm,lv in allTreeLeafs.items() if lvNm.startswith(prefix) and lvNm != sizeNm and lv.GetLeafCount().GetName() == sizeNm)
        sizeOp = GetColumn(allTreeLeafs[sizeNm].GetTypeName(), sizeNm)
        if allTreeLeafs[sizeNm].GetTypeName() != SizeType:
            sizeOp = adaptArg(op.static_cast(SizeType, sizeOp))
        for lvNm in itm_lvs:
            lvNm_short = lvNm[len(prefix):]
            col = GetArrayColumn(allTreeLeafs[lvNm].GetTypeName(), lvNm, sizeOp).result
            itm_dict[lvNm_short] = itemProxy(col)
        ## create p4 branches
        p4AttNames = ("pt", "eta", "phi", "mass")
        if all(("".join((prefix, att)) in itm_lvs) for att in p4AttNames):
            itm_dict["p4"] = funProxy(lambda inst : Construct("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >", (inst.pt, inst.eta, inst.phi, inst.mass)).result)
            ## mother/daughter references
        itm_bases = [ContainerGroupItemProxy]
        if grpNm == "genpart":
            itm_bases.append(Phase2GenRanges)
        itmcls = type("{0}GroupItemProxy".format(grpNm.capitalize()), tuple(itm_bases), itm_dict)
        coll = ContainerGroupProxy(prefix, None, sizeOp, itmcls)
        tree_dict[grpNm] = coll
        tree_children.append(coll)

        for lvNm in itm_lvs:
            del tree_dict[lvNm]
        del tree_dict[sizeNm] ## go through op.rng_len

    TreeProxy = type("{0}Proxy".format(aTree.GetName().capitalize()), (TreeBaseProxy,), tree_dict)
    treeProxy = TreeProxy(aTree)

    for pc in tree_children:
        pc._parent = treeProxy

    return treeProxy
