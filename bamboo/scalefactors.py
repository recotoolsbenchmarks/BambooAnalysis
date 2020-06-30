"""
The :py:mod:`bamboo.scalefactors` module contains helper methods
for configuring scale factors, fake rates etc.

The basic configuration parameter is the JSON file path for a set of scalefactors.
There two basic types are

- lepton scale factors (dependent on a number of object variables, e.g. pt and eta),
- jet (b-tagging) scale factors (grouped set for different flavours, for convenience)

Different values (depending on the data-taking period)
can be taken into account by weighting or by randomly sampling.
"""
__all__ = ("get_scalefactor", "lumiPerPeriod")

from itertools import chain
from functools import partial

from . import treefunctions as op

#: Integrated luminosity (in 1/pb) per data taking period
lumiPerPeriod_default = {
    ## 2016 - using approved normtag + 07Aug2017 re-reco golden JSON
    # ReReco/Final/Cert_271036-284044_13TeV_ReReco_07Aug2017_Collisions16_JSON.txt
      "Run2016B" : 5750.491
    , "Run2016C" : 2572.903
    , "Run2016D" : 4242.292
    , "Run2016E" : 4025.228
    , "Run2016F" : 3104.509
    , "Run2016G" : 7575.824
    , "Run2016H" : 8650.628
    # hww muon periods
    , "Run271036to275783" : 6274.191
    , "Run275784to276500" : 3426.131
    , "Run276501to276811" : 3191.207

    ## 2017 - using approved normtag + re-reco golden JSON
    # ReReco/Cert_294927-306462_13TeV_EOY2017ReReco_Collisions17_JSON_v1.txt
    , "Run2017B" : 4793.970
    , "Run2017C" : 9632.746
    , "Run2017D" : 4247.793
    , "Run2017E" : 9314.581
    , "Run2017F" : 13539.905
    
    ## 2018 - using approved normtag + re-reco golden JSON
    # ReReco/Cert_314472-325175_13TeV_17SeptEarlyReReco2018ABC_PromptEraD_Collisions18_JSON.txt
    , "Run2018A" : 14027.614
    , "Run2018B" : 7066.552
    , "Run2018C" : 6898.817
    , "Run2018D" : 31747.582
    # before/after muon HLT update (during 2018A)
    , "Run315264to316360" : 8928.970
    , "Run316361to325175" : 50789.746
    }


# TODO maybe move this elsewhere
binningVariables_nano = {
      "Eta" : lambda obj : obj.eta
    , "AbsEta" : lambda obj : op.abs(obj.eta)
    , "ClusEta" : lambda el : el.eta+el.deltaEtaSC
    , "AbsClusEta" : lambda el : op.abs(el.eta+el.deltaEtaSC)
    , "Pt" : lambda obj : obj.pt
    }

def getBinningVarNames(jsonpath):
    import json
    with open(jsonpath, "r") as jsf:
        cont = json.load(jsf)
    return tuple(cont["variables"])

class BinningParameters(object):
    def __init__(self, binningVars):
        self.binningVars = binningVars
    def __call__(self, obj):
        return op.construct("Parameters",
                       (op.initList("std::initializer_list<Parameters::value_type::value_type>", "Parameters::value_type::value_type", (
                           op.initList("Parameters::value_type::value_type", "float", (op.extVar("int", "BinningVariable::{0}".format(bvNm.replace("ClusEta", "Eta"))), bv(obj)))
                           for bvNm,bv in self.binningVars.items())),)
                   )

def getBinningParameters(bVarNames, isElectron=False, moreVars=dict(), paramDefs=dict()):
    if isElectron:
        bVarNames = [ k.replace("Eta", "ClusEta") for k in bVarNames ]
    theDict = dict(paramDefs)
    theDict.update(moreVars)
    return BinningParameters(dict((k,theDict[k]) for k in bVarNames))

class ScaleFactor(object):
    def __init__(self, cppDef=None, args=None, iface="ILeptonScaleFactor", systName=None, seedFun=None):
        self._cppDef = cppDef
        self._args = args
        self.sfOp = op.define(iface, cppDef)
        self._systName = systName
        self._seedFun = seedFun
    def __call__(self, obj, variation="Nominal"):
        from .treedecorators import makeConst, boolType
        from .treeoperations import ScaleFactorWithSystOp
        expr = self.sfOp.get(*tuple(chain(
                   list(a(obj) for a in self._args)
                 , ([ self._seedFun(obj) ] if self._seedFun else [])
                 , (op.extVar("int", variation),)
               )))
        if self._systName and variation == "Nominal": ## wrap
            expr._parent = ScaleFactorWithSystOp(expr._parent, self._systName)
        return expr

def get_scalefactor(objType, key, combine=None, additionalVariables=dict(), sfLib=dict(), paramDefs=dict(), lumiPerPeriod=dict(), periods=None, getFlavour=None, isElectron=False, systName=None, seedFun=None):
    """ Construct a scalefactor callable

    :param objType: object type: ``"lepton"``, ``"dilepton"``, or ``"jet"``
    :param key: key in ``sfLib`` (or tuple of keys, in case of a nested dictionary)
    :param sfLib: dictionary (or nested dictionary) of scale factors. A scale factor entry is either a path to a JSON file, or a list of pairs ``(periods, path)``, where ``periods`` is a list of periods found in ``lumiPerPeriod`` and ``path`` is a path to the JSON file with the scale factors corresponding to those run periods.
    :param combine: combination strategy for combining different run periods (``"weight"`` or ``"sample"``)
    :param paramDefs: dictionary of binning variable definitions (name to callable)
    :param additionalVariables: additional binning variable definitions (TODO: remove)
    :param lumiPerPeriod: alternative definitions and relative weights of run periods
    :param periods: Only combine scale factors for those periods
    :param isElectron: if True, will use supercluster eta instead of eta (for ``"lepton"`` type only) (TODO: find a better way of doing that)
    :param systName: name of the associated systematic nuisance parameter
    :param seedFun: (only when combining scalefactor by sampling) callable to get a random generator seed for an object, e.g. ``lambda l : l.idx+42``

    :returns: a callable that takes ``(object, variation="Nominal")`` and returns a floating-point number proxy
    """
    ##
    ## Interpret args, get defaults etc.
    ## 
    if isinstance(key, tuple):
        # interpret key=("a", "b") as ...["a"]["b"]
        mainKey = key[0]
        config = sfLib[key[0]]
        for idx in range(1,len(key)):
            config = config[key[idx]]
    else:
        mainKey = key
        config = sfLib[key]

    if combine is not None:
        combPrefix = { "weight" : "W"
                     , "sample" : "Smp" }.get(combine, "W")
        if combine == "sample" and not seedFun:
            raise ValueError("If combining by sampling, a seed function needs to be passed to get_scalefactor")

    if getFlavour is None:
        getFlavour = lambda j : j.hadronFlavour
    getFlavour = partial(lambda getter, j : op.extMethod("IJetScaleFactor::get_flavour")(getter(j)), getFlavour)

    lumiPerPeriod_default.update(lumiPerPeriod)
    lumiPerPeriod = lumiPerPeriod_default
    if periods is None:
        periods = lumiPerPeriod.keys()

    ##
    ## Construct scalefactors
    ##
    if objType == "lepton":
        iface = "ILeptonScaleFactor"
        if isinstance(config, str):
            return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(config),
                    args=(getBinningParameters(getBinningVarNames(config), isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                    iface=iface, systName=systName)
        else:
            if combPrefix == "":
                raise ValueError("A combination mode needs to be specified for this scale factor")
            selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),path)
                    for eras,path in config if any(ier in periods for ier in eras))))
            if len(selConfigs) < 1:
                raise RuntimeError("Zero period configs selected for config {0} with periods {1}".format(", ".join("({0} : {1})".format(list(eras), path) for eras, path in config), list(periods)))
            elif len(selConfigs) == 1:
                return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(selConfigs[0][1]),
                        args=(getBinningParameters(getBinningVarNames(selConfigs[0][1]), isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                        iface=iface, systName=systName)
            else:
                bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,iPth in selConfigs))
                return ScaleFactor(cppDef=(
                            'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                'std::make_unique<ScaleFactor>("{0}")'.format(path) for wgt, path in selConfigs), iface=iface)+
                            'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,path in selConfigs), cmb=combPrefix)+
                              'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                            ),
                        args=(getBinningParameters(bVarNames, isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                        iface=iface, systName=systName, seedFun=(seedFun if combine == "sample" else None))
    elif objType == "dilepton":
        iface = "IDiLeptonScaleFactor"
        if isinstance(config, tuple) and len(config) == 4:
            if not all(isinstance(iCfg, str) for iCfg in config):
                raise TypeError("Config for dilepton scale factor should be quadruplet of paths or list f weights and triplets, found {0}".format(config))

            return ScaleFactor(cppDef="const DiLeptonFromLegsScaleFactor <<name>>{{{0}}};".format(", ".join(
                        'std::make_unique<ScaleFactor>("{0}")'.format(leplepCfg) for leplepCfg in config)),
                    args=[ (lambda bp : (lambda ll : bp(ll[0])))(getBinningParameters(set(chain(getBinningVarNames(config[0]), getBinningVarNames(config[1]))), moreVars=additionalVariables, paramDefs=paramDefs))
                         , (lambda bp : (lambda ll : bp(ll[1])))(getBinningParameters(set(chain(getBinningVarNames(config[2]), getBinningVarNames(config[3]))), moreVars=additionalVariables, paramDefs=paramDefs)) ],
                    iface=iface, systName=systName)
        else:
            raise NotImplementedError("Still to do this part")
    elif objType == "jet":
        iface = "IJetScaleFactor"
        if isinstance(config, tuple) and len(config) == 3:
            if not all(isinstance(iCfg, str) for iCfg in config):
                raise TypeError("Config for b-tagging should be triplet of paths or list of weights and triplets, found {0}".format(config))
            else:
                bVarNames = set(chain.from_iterable(getBinningVarNames(iCfg) for iCfg in config))
                return ScaleFactor(cppDef='const BTaggingScaleFactor <<name>>{{{0}}};'.format(", ".join('"{0}"'.format(iCfg) for iCfg in config)),
                        args=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), getFlavour),
                        iface=iface, systName=systName)
        else:
            if not ( all((isinstance(iCfg, tuple) and len(iCfg) == 3 and all(isinstance(iPth, str) for iPth in iCfg) ) for iCfg in config) ):
                raise TypeError("Config for b-tagging should be triplet of paths or list of weights and triplets, found {0}".format(config))
            else:
                if combPrefix == "":
                    raise ValueError("A combination mode needs to be specified for this scale factor")
                selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                    ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),path)
                        for eras,path in config if any(ier in periods for ier in eras))))
                if len(selConfigs) < 1:
                    raise RuntimeError("Zero configs")
                elif len(selConfigs) == 1:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iCfg) for iCfg in selConfigs[0]))
                    return ScaleFactor(cppDef='const BTaggingScaleFactor <<name>>{{{0}}};'.format(", ".join('"{0}"'.format(iCfg) for iCfg in selConfigs[0])),
                            args=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), getFlavour),
                            iface=iface, systName=systName)
                else:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,paths in selConfigs for iPth in paths))
                    return ScaleFactor(cppDef=(
                                'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                    'std::make_unique<BTaggingScaleFactor>({0})'.format(", ".join('"{0}"'.format(iPth) for iPth in paths)) for wgt, paths in selConfigs), iface=iface)+
                                'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,paths in selConfigs), cmb=combPrefix)+
                                  'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                                ),
                            arg=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), getFlavour),
                            iface=iface, systName=systName, seedFun=(seedFun if combine == "sample" else None))
    else:
        raise ValueError("Unknown object type: {0}".format(objType))

class BtagSF:
    def _nano_getPt(jet):
        import bamboo.treeproxies as _tp
        if isinstance(jet._parent, _tp.AltCollectionProxy):
            bs = jet._parent._base
            return bs.brMap["pt"].wrapped.result[jet.idx] ## use nominal always
        else:
            return jet.pt
    def _nano_getEta(jet):
        import bamboo.treefunctions as op
        return op.abs(jet.eta)
    def _nano_getJetFlavour(jet):
        import bamboo.treefunctions as op
        return op.extMethod("BTagEntry::jetFlavourFromHadronFlavour")(jet.hadronFlavour)

    def _translate_btagSFVarToJECVar(btagVarName):
        if btagVarName.startswith("up_jes") or btagVarName.startswith("down_jes"):
            return f"jes{1}{0}".format(*btagVarName.split("_jes"))
        else:
            return btagVarName

    def __init__(self, taggerName, csvFileName, wp=None, sysType="central", otherSysTypes=None, systName="jet", measurementType=None, getters=None, jesTranslate=None, sel=None, uName=None):
        """
        Declare a BTagCalibration (if needed) and BTagCalibrationReader (unique, based on ``uName``), and decorate for evaluation

        :param taggerName: first argument for ``BTagCalibration``
        :param csvFileName: name of the CSV file with scalefactors
        :param wp: working point (used as ``BTagEntry::OP_{wp.upper()}``)
        :param sysType: nominal value systematic type (``"central"``, by default)
        :param otherSysTypes: other systematic types to load in the reader
        :param systName: name of the systematic ("jet", by default, to coordinate the JEC variations)
        :param measurementType: dictionary with measurement type per true flavour, or a string if the same for all (if not specified, ``"comb"`` will be used for b- and c-jets, and ``incl`` for light-flavour jets)
        :param getters: dictionary of methods to get the kinematics and classifier for a jet (the keys ``Pt``, ``Eta``, ``JetFlavour``, and ``Discri`` are used. For the former three, the defaults are for NanoAOD)
        :param jesTranslate: translation function for JEC systematic variations, from the names in the CSV file to those used for the jets (the default default should work for on-the-fly corrections)
        :param sel: a selection in the current graph
        :param uName: unique name, to declare the reader (e.g. sample name)
        """
        if otherSysTypes is None:
            otherSysTypes = []
        self.nomName = sysType
        self.varNames = otherSysTypes
        self.systName = systName
        if measurementType is None: ## BTV recommendation for b-tagging with fixed working points
            measurementType = {"B": "comb", "C": "comb", "UDSG": "incl"}
        elif isinstance(measurementType, str): ## if string, use it for all
            measurementType = {fl: measurementType for fl in ("B", "C", "UDSG")}
        calibName = sel._fbe.symbol(f'const BTagCalibration <<name>>{{"{taggerName}", "{csvFileName}"}};', nameHint=f"bTagCalib_{taggerName}")
        readerName = sel._fbe.symbol('BTagCalibrationReader <<name>>{{BTagEntry::OP_{0}, "{1}", {{ {2} }} }}; // for {3}'.format(wp.upper(), sysType, ", ".join(f'"{sv}"' for sv in otherSysTypes), uName), nameHint="bTagReader_{0}".format("".join(c for c in uName if c.isalnum())))
        from bamboo.root import gbl
        calibHandle = getattr(gbl, calibName)
        readerHandle = getattr(gbl, readerName)
        for flav,measType in measurementType.items():
            readerHandle.load(calibHandle, getattr(gbl.BTagEntry, f"FLAV_{flav}"), measType)
        import bamboo.treefunctions as op
        self.reader = op.extVar("BTagCalibrationReader", readerName)
        if getters is None:
            self.getters = [ BtagSF._nano_getJetFlavour, BtagSF._nano_getEta, BtagSF._nano_getPt ]
        else:
            self.getters = [ getters.get("JetFlavour", BtagSF._nano_getJetFlavour), getters.get("Eta", BtagSF._nano_getEta), getters.get("Pt", BtagSF._nano_getPt) ]
            if "Discri" in getters:
                self.getters.append(getters["Discri"])
        self.jesTranslate = jesTranslate if jesTranslate is not None else BtagSF._translate_btagSFVarToJECVar
    def _evalFor(self, var, jet):
        import bamboo.treefunctions as op
        return self.reader.eval_auto_bounds(op._tp.makeConst(var, "std::string"), *(gtr(jet) for gtr in self.getters))
    def __call__(self, jet, systName=None, nomVar=None, systVars=None):
        """
        Evaluate the scalefactor for a jet

        Please note that this only gives the applicable scalefactor:
        to obtain the event weight one of the recipes in the
        `POG twiki<https://twiki.cern.ch/twiki/bin/viewauth/CMS/BTagSFMethods>`_
        should be used.

        By default the systematic name, nominal and systematic variations
        are taken from the :py:class:`bamboo.scalefactors.BtagSF` instances,
        but they can be overriden with the ``systName``, ``nomVar`` and ``systVars``
        keyword arguments.
        """
        import bamboo.treefunctions as op
        nom = self._evalFor((nomVar if nomVar is not None else self.nomName), jet)
        if systVars is not None and not systVars: ## empty, no systematic
            return nom
        else:
            if systVars is not None:
                if any(sv not in self.varNames for sv in systVars):
                    logger.error("B-tag SF systematic variations {0} were not loaded, skipping".format(", ".join(sv for sv in systVars if sv not in self.varNames)))
                    systVars = [ sv for sv in systVars if sv in self.varNames ]
            else:
                systVars = self.varNames
            return op.systematic(nom, name=(systName if systName is not None else self.systName),
                    **{ self.jesTranslate(var) : self._evalFor(var, jet) for var in systVars })
