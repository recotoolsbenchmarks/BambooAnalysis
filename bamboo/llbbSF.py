"""
2016 scalefactors from the cp3-llbb HH and ZA analyses
"""
from . import treefunctions as op
from . import scalefactors

from itertools import chain

def localize_llbbFwk(aPath, package="Framework", hat="ScaleFactors"):
    import os.path
    if "CP3LLBBBASE" not in os.environ:
        raise RuntimeError("This version of the scalefactors needs $CP3LLBBASE to find the scalefactor files (to be improved)")
    return os.path.join(os.environ["CP3LLBBBASE"], package, "data", hat, aPath)

binningVariables = {
      "Eta"       : lambda obj : obj.p4.Eta()
    , "ClusEta"   : lambda obj : obj.clusterEta
    , "AbsEta"    : lambda obj : op.abs(obj.p4.Eta())
    , "AbsClusEta": lambda obj : op.abs(obj.clusterEta)
    , "Pt"        : lambda obj : obj.p4.Pt()
    } ## TODO add more?

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

def get_scalefactor(objType, key, combine=None, periods=None, isElectron=False, additionalVariables=dict()):
    return scalefactors.get_scalefactor(objType, key, combine=combine, periods=periods, additionalVariables=dict(), sfLib=all_scalefactors, paramDefs=binningVariables, isElectron=isElectron)

if __name__ == "__main__":
    lSF = get_scalefactors("lepton", ("muon_2016_80", "iso_loose_id_loose"), combine="weight")
    cmva_discriVar = {"BTagDiscri":lambda j : j.pfCombinedMVAV2BJetTags}
    jSF = get_scalefactors("jet", ("btag_2016_moriond2017", "cMVAv2_loose"), combine="weight", additionalVariables=cmva_discriVar)
    #llSF = get_scalefactors("dilepton", "eleltrig_2016_HHMoriond17")
