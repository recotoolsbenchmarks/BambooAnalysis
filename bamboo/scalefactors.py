"""
Helpers for configuring scale factors, fake rates etc.

The basic configuration parameter is the json file path for a set of factors.
There two basic types are
- lepton scale factors (dependent on a number of object variables, e.g. PT and ETA),
- jet (b-tagging) scale factors (grouped set for different flavours, for convenience)

Different values (depending on the data-taking period)
can be taken into account by weighting or by randomly sampling.
"""
__all__ = ("get_scalefactor",)

from itertools import chain

def localize_llbbFwk(aPath, package="Framework", hat="ScaleFactors"):
    import os.path
    if "CP3LLBBBASE" not in os.environ:
        raise RuntimeError("This version of the scalefactors needs $CP3LLBBASE to find the scalefactor files (to be improved)")
    return os.path.join(os.environ["CP3LLBBBASE"], package, "data", hat, aPath)


lumiPerPeriod = {
      "Run2016B" : 5785.152 ## averaged 5783.740 (DoubleMuon), 5787.976 (DoubleEG) and 5783.740 (MuonEG) - max dev. from average is 5.e-4 << lumi syst
    , "Run2016C" : 2573.399
    , "Run2016D" : 4248.384
    , "Run2016E" : 4009.132
    , "Run2016F" : 3101.618
    , "Run2016G" : 7540.488
    , "Run2016H" : 8605.689
    ##
    , "Run271036to275783" : 6274.191
    , "Run275784to276500" : 3426.131
    , "Run276501to276811" : 3191.207
    }

hwwMuonPeriods_2016 = [ "Run271036to275783", "Run275784to276500", "Run276501to276811" ]
hwwElePeriods_2016 = [] ## TODO

all_scalefactors = {
      "electron_2015_76"  : dict((k,localize_llbbFwk(v)) for k, v in chain(
          dict(("id_{wp}".format(wp=wp.lower()), ("Electron_CutBasedID_{wp}WP_fromTemplates_withSyst_76X.json".format(wp=wp)))
         for wp in ("Veto", "Loose", "Medium", "Tight")).items()
        , { "hww_wp" : "Electrons_HWW_CutBasedID_TightWP_76X_forHWW_Final.json" }.items()
        ))
    , "muon_2015_76" : dict((k, localize_llbbFwk(v)) for k, v in chain(
          dict(("id_{wp}".format(wp=wp.lower()), "Muon_{wp}ID_genTracks_id.json".format(wp=wp)) for wp in ("Soft", "Loose", "Medium")).items()
        , { "id_tight" : "Muon_TightIDandIPCut_genTracks_id.json" }.items()
        , dict(("iso_{isowp}_id_{idwp}".format(isowp=isowp.lower(), idwp=idwp.lower()), "Muon_{isowp}RelIso_{idwp}ID_iso.json".format(isowp=isowp, idwp=idwp))
            for (idwp,isowp) in (("Loose", "Loose"), ("Loose", "Medium"), ("Loose", "Tight"), ("Tight", "Medium"), ("Tight", "Tight")) ).items()
        , { "id_hww"   : "Muon_MediumID_Data_MC_25ns_PTvsETA_HWW_76.json"
          , "iso_tight_id_hww" : "Muon_ISOTight_Data_MC_25ns_PTvsETA_HWW.json" }.items()
        ))
    , "btag_2015_76" : dict() ## TODO
    ## 2016 CMSSW_8_0_...
    , "muon_2016_80" : dict((k,( localize_llbbFwk(v) if isinstance(v, str)
                               else [ (eras, localize_llbbFwk(path)) for eras,path in v ]))
                           for k, v in chain(
        ## https://twiki.cern.ch/twiki/bin/view/CMS/MuonReferenceEffsRun2#Tracking_efficiency_provided_by
          { "tracking" : "Muon_tracking_BCDEFGH.json" }.items()
        ## https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonWorkInProgressAndPagResults
        , dict(("id_{wp}".format(wp=wp.lower()), [ (("Run2016{0}".format(ltr) for ltr in eras), "Muon_{wp}ID_genTracks_id_{era}.json".format(wp=wp, era=eras)) for eras in ("BCDEF", "GH") ]) for wp in ("Loose", "Medium", "Tight")).items()
        , { "id_medium2016" : [ (("Run2016{0}".format(ltr) for ltr in eras), "Muon_MediumID2016_genTracks_id_{era}.json".format(era=eras)) for eras in ("BCDEF", "GH") ] }.items()
        , dict(("iso_{isowp}_id_{idwp}".format(isowp=isowp.lower(), idwp=idwp.lower())
               , [ (("Run2016{0}".format(ltr) for ltr in eras), "Muon_{isowp}ISO_{idwp}ID_iso_{era}.json".format(isowp=isowp, idwp=idwp, era=eras)) for eras in ("BCDEF", "GH") ])
            for (isowp, idwp) in [("Loose", "Loose"), ("Loose", "Medium"), ("Loose", "Tight"), ("Tight", "Medium"), ("Tight", "Tight")]).items()
        ## https://twiki.cern.ch/twiki/bin/view/CMS/HWW2016TriggerAndIdIsoScaleFactorsResults
        , { "iso_tight_hww" : [ ((era,), "Muon_data_mc_ISOTight_Run2016_{era}_PTvsETA_HWW.json".format(era=era)) for era in hwwMuonPeriods_2016 ] }.items()
        , { "id_tight_hww"  : [ ((era,), "Muon_data_mc_TightID_Run2016_{era}_PTvsETA_HWW.json".format(era=era)) for era in hwwMuonPeriods_2016 ] }.items()
        ))
    ## https://twiki.cern.ch/twiki/bin/view/CMS/EgammaIDRecipesRun2#Efficiencies_and_scale_factors
    , "electron_2016_ichep2016" : dict((k,( localize_llbbFwk(v) if isinstance(v, str)
                               else [ (eras, localize_llbbFwk(path)) for eras,path in v ]))
                           for k, v in chain(
          { "gsf_tracking"  : "Electron_EGamma_SF2D_gsfTracking.json" }.items()
        , dict(("id_{wp}".format(wp=wp), "Electron_EGamma_SF2D_{wp}.json".format(wp=wp)) for wp in ("veto", "loose", "medium", "tight")).items()
        , { "hww_wp"        : [ ((era,), "Electron_Tight_Run2016_{era}_HWW.json".format(era=era)) for era in hwwElePeriods_2016 ] }.items()
        ))
    , "electron_2016_moriond2017" : dict((k,( localize_llbbFwk(v) if isinstance(v, str)
                               else [ (eras, localize_llbbFwk(path)) for eras,path in v ]))
                           for k, v in chain(
          { "gsf_tracking"  : "Electron_EGamma_SF2D_reco_moriond17.json" }.items()
        , dict(("id_{wp}".format(wp=wp), "Electron_EGamma_SF2D_{wp}_moriond17.json".format(wp=wp)) for wp in ("veto", "loose", "medium", "tight")).items()
        ))
    , "eleltrig_2016_HHMoriond17" : tuple(localize_llbbFwk("{0}.json".format(nm), package="ZAAnalysis", hat="Efficiencies") for nm in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"))
    , "elmutrig_2016_HHMoriond17" : tuple(localize_llbbFwk("{0}.json".format(nm), package="ZAAnalysis", hat="Efficiencies") for nm in ("Electron_IsoEle23Leg", "Electron_IsoEle12Leg", "Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg"))
    , "mueltrig_2016_HHMoriond17" : tuple(localize_llbbFwk("{0}.json".format(nm), package="ZAAnalysis", hat="Efficiencies") for nm in ("Muon_XPathIsoMu23leg", "Muon_XPathIsoMu8leg", "Electron_IsoEle23Leg", "Electron_IsoEle12Leg"))
    , "mumutrig_2016_HHMoriond17" : tuple(localize_llbbFwk("{0}.json".format(nm), package="ZAAnalysis", hat="Efficiencies") for nm in ("Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg", "Muon_DoubleIsoMu17Mu8_IsoMu17leg", "Muon_DoubleIsoMu17TkMu8_IsoMu8legORTkMu8leg"))
    ## https://twiki.cern.ch/twiki/bin/view/CMS/BtagRecommendation80XReReco
    , "btag_2016_moriond2017" : dict((k,( tuple(localize_llbbFwk(fv) for fv in v) if isinstance(v,tuple) and all(isinstance(fv, str) for fv in v)
                               else [ (eras, tuple(localize_llbbFwk(fpath) for fpath in paths)) for eras,paths in v ]))
                           for k, v in 
          dict(("cMVAv2_{wp}".format(wp=wp), tuple("BTagging_{wp}_{flav}_{calib}_cmvav2_BtoH_moriond17.json".format(wp=wp, flav=flav, calib=calib) for (flav, calib) in (("lightjets", "incl"), ("cjets", "ttbar"), ("bjets", "ttbar")))) for wp in ("loose", "medium", "tight")).items()
        )
    }

from . import treefunctions as op
binningVariablesByName = {
      "Eta"       : lambda obj : obj.p4.Eta()
    , "ClusEta"   : lambda obj : obj.clusterEta
    , "AbsEta"    : lambda obj : op.abs(obj.p4.Eta())
    , "AbsClusEta": lambda obj : op.abs(obj.clusterEta)
    , "Pt"        : lambda obj : obj.p4.Pt()
    } ## TODO add more?

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

def getBinningParameters(bVarNames, isElectron=False, moreVars=dict()):
    if isElectron:
        bVarNames = [ k.replace("Eta", "ClusEta") for k in bVarNames ]
    theDict = dict(binningVariablesByName)
    theDict.update(moreVars)
    return BinningParameters(dict((k,theDict[k]) for k in bVarNames))

class ScaleFactor(object):
    def __init__(self, cppDef=None, args=None, iface="ILeptonScaleFactor"):
        self._cppDef = cppDef
        self._args = args
        self.sfOp = op.define(iface, cppDef)
    def __call__(self, obj, variation="Nominal", withMCCheck=True):
        from .treedecorators import makeConst, boolType
        expr = self.sfOp.get(*tuple(chain(
                   list(a(obj) for a in self._args)
                 , (op.extVar("int", variation),)
               )))
        if withMCCheck:
            expr = op.switch(op.extVar(boolType, "runOnMC"), expr, makeConst(1., "float")) ## TODO declare runOnMC then
        return expr

def get_scalefactor(objType, key, periods=None, combine=None, additionalVariables=dict()):
    ##
    ## Interpret args, get defaults etc.
    ## 
    if isinstance(key, tuple):
        # interpret key=("a", "b") as ...["a"]["b"]
        mainKey = key[0]
        config = all_scalefactors[key[0]]
        for idx in range(1,len(key)):
            config = config[key[idx]]
    else:
        mainKey = key
        config = all_scalefactors[key]

    if periods is None:
        if "2016" in mainKey:
            periods = [ "Run2016{0}".format(ltr) for ltr in "BCDEFGH" ]
        else:
            periods = []
    periods = set(periods)
    
    if combine is not None:
        combPrefix = { "weight" : "W"
                     , "sample" : "Smp" }.get(combine, "W")

    ##
    ## Construct scalefactors
    ##
    if objType == "lepton":
        iface = "ILeptonScaleFactor"
        isElectron = (key[0].split("_")[0] == "electron")
        if isinstance(config, str):
            return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(config),
                    args=(getBinningParameters(getBinningVarNames(config), isElectron=isElectron, moreVars=additionalVariables),),
                    iface=iface)
        else:
            if combPrefix == "":
                raise ValueError("A combination mode needs to be specified for this scale factor")
            selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),path)
                    for eras,path in config if any(ier in periods for ier in eras))))
            if len(selConfigs) < 1:
                raise RuntimeError("Zero configs")
            elif len(selConfigs) == 1:
                return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(selConfigs[0][1]),
                        args=(getBinningParameters(getBinningVarNames(selConfigs[0][1]), isElectron=isElectron, moreVars=additionalVariables),),
                        iface=iface)
            else:
                bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,iPth in selConfigs))
                return ScaleFactor(cppDef=(
                            'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                'std::make_unique<ScaleFactor>("{0}")'.format(path) for wgt, path in selConfigs), iface=iface)+
                            'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,path in selConfigs), cmb=combPrefix)+
                              'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                            ),
                        args=(getBinningParameters(bVarNames, isElectron=isElectron, moreVars=additionalVariables),),
                        iface=iface)
    elif objType == "dilepton":
        iface = "IDiLeptonScaleFactor"
        if isinstance(config, tuple) and len(config) == 4:
            if not all(isinstance(iCfg, str) for iCfg in config):
                raise TypeError("Config for dilepton scale factor should be quadruplet of paths or list f weights and triplets, found {0}".format(config))

            return ScaleFactor(cppDef="const DiLeptonFromLegsScaleFactor <<name>>{{{0}}};".format(", ".join(
                        'std::make_unique<LeptonScaleFactor>("{0}")'.format(leplepCfg) for leplepCfg in config)),
                    args=[ (lambda bp : (lambda ll : bp(ll[0])))(getBinningParameters(set(chain(getBinningVarNames(config[0]), getBinningVarNames(config[1]))), moreVars=additionalVariables))
                         , (lambda bp : (lambda ll : bp(ll[1])))(getBinningParameters(set(chain(getBinningVarNames(config[2]), getBinningVarNames(config[3]))), moreVars=additionalVariables)) ],
                    iface=iface)
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
                        args=(getBinningParameters(bVarNames, moreVars=additionalVariables), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                        iface=iface)
        else:
            if not ( all((isinstance(iCfg, tuple) and len(iCfg) == 3 and all(isinstance(iPth, str) for iPth in iCfg) ) for iCfg in config) ):
                raise TypeError("Config for b-tagging should be triplet of paths or list of weights and triplets, found {0}".format(config))
            else:
                if combPrefix == "":
                    raise ValueError("A combination mode needs to be specified for this scale factor")
                selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                    ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),paths)
                        for eras,paths in config if any(ier in periods for ier in eras))))
                if len(selConfigs) < 1:
                    raise RuntimeError("Zero configs")
                elif len(selConfigs) == 1:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iCfg) for iCfg in selConfigs[0]))
                    return ScaleFactor(cppDef='const BTaggingScaleFactor <<name>>{{{0}}};'.format(", ".join('"{0}"'.format(iCfg) for iCfg in selConfigs[0])),
                            args=(getBinningParameters(bVarNames, moreVars=additionalVariables), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                            iface=iface)
                else:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,paths in selConfigs for iPth in paths))
                    return ScaleFactor(cppDef=(
                                'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                    'std::make_unique<BTaggingScaleFactor>({0})'.format(", ".join('"{0}"'.format(iPth) for iPth in paths)) for wgt, paths in selConfigs), iface=iface)+
                                'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,paths in selConfigs), cmb=combPrefix)+
                                  'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                                ),
                            arg=(getBinningParameters(bVarNames, moreVars=additionalVariables), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                            iface=iface)
    else:
        raise ValueError("Unknown object type: {0}".format(objType))

if __name__ == "__main__":
    lSF = get_scalefactors("lepton", ("muon_2016_80", "iso_loose_id_loose"), combine="weight")
    cmva_discriVar = {"BTagDiscri":lambda j : j.pfCombinedMVAV2BJetTags}
    jSF = get_scalefactors("jet", ("btag_2016_moriond2017", "cMVAv2_loose"), combine="weight", additionalVariables=cmva_discriVar)
    #llSF = get_scalefactors("dilepton", "eleltrig_2016_HHMoriond17")
