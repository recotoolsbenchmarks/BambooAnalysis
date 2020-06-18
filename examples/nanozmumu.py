"""
Example analysis module: make a dimuon mass plot from a NanoAOD
"""
from bamboo.analysismodules import NanoAODModule, NanoAODHistoModule, NanoAODSkimmerModule
import os.path
from bamboo.logging import getLogger
logger = getLogger(__name__)

import bamboo.plots
class Plot(bamboo.plots.Plot):
    def produceResults(self, bareResults, fbe):
        if any("__qcdScale" in h.GetName() for h in bareResults):
            hNom = next(h for h in bareResults if "__" not in h.GetName())
            prefix = f"{hNom.GetName()}__qcdScale"
            hVar_qcdScale = [ h for h in bareResults if h.GetName().startswith(prefix) ]
            if not all(hv.GetNcells() == hNom.GetNcells() for hv in hVar_qcdScale):
                logger.error("Variation histograms do not have the same binning as the nominal histogram")
            elif len(hVar_qcdScale) < 2:
                logger.error("At least two variations histograms must be provided")
            else: ## make an envelope from maximum deviations
                import numpy as np
                vars_cont = np.array([ [ hv.GetBinContent(i) for i in range(hv.GetNcells()) ] for hv in hVar_qcdScale ])
                hVar_up = hNom.Clone(f"{prefix}up")
                hVar_down = hNom.Clone(f"{prefix}down")
                from itertools import count
                for i,vl,vh in zip(count(), np.amin(vars_cont, axis=0), np.amax(vars_cont, axis=0)):
                    hVar_down.SetBinContent(i, vl)
                    hVar_up.SetBinContent(i, vh)
                return bareResults + [ hVar_up, hVar_down ]
        return bareResults

def makeBTagCalibrationReader(taggerName, csvFileName, wp=None, sysType="central", otherSysTypes=None, measurementType=None, sel=None, uName=None):
    """
    Declare a BTagCalibration (if needed) and BTagCalibrationReader (unique, based on ``uName``)

    :param taggerName: first argument for ``BTagCalibration``
    :param csvFileName: name of the CSV file with scalefactors
    :param wp: working point (used as ``BTagEntry::OP_{wp.upper()}``)
    :param systType: nominal value systematic type (``"central"``, by default)
    :param otherSysTypes: other systematic types to load in the reader
    :param measurementType: dictionary with measurement type per true flavour, or a string if the same for all (if not specified, ``"comb"`` will be used for b- and c-jets, and ``incl`` for light-flavour jets)
    :param sel: a selection in the current graph
    :param uName: unique name, to declare the reader (e.g. sample name)
    """
    if otherSysTypes is None:
        otherSysTypes = []
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
    return op.extVar("BTagCalibrationReader", readerName)

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

    def __init__(self, reader, getPt=None, getEta=None, getJetFlavour=None, getDiscri=None):
        self.reader = reader
        self.getPt = getPt if getPt is not None else BtagSF._nano_getPt
        self.getEta = getEta if getEta is not None else BtagSF._nano_getEta
        self.getJetFlavour = getJetFlavour if getJetFlavour is not None else BtagSF._nano_getJetFlavour
        self.getDiscri = getDiscri

    def _evalFor(self, var, jet):
        import bamboo.treefunctions as op
        args = [ op._tp.makeConst(var, "std::string"), self.getJetFlavour(jet), self.getEta(jet), self.getPt(jet) ]
        if self.getDiscri:
            args.append(self.getDiscri(jet))
        return self.reader.eval_auto_bounds(*args)

    def __call__(self, jet, systVars=None, systName=None):
        import bamboo.treefunctions as op
        nom = self._evalFor("central", jet)
        if systVars is None:
            return nom
        else:
            if systName is None:
                raise RuntimeError("A name for the systematic uncertainty group is needed")
            return op.systematic(nom, name=systName, **{ var : self._evalFor(var, jet) for var in systVars })

class NanoZMuMuBase(NanoAODModule):
    """ Base module for NanoAOD Z->MuMu example """
    def addArgs(self, parser):
        super(NanoZMuMuBase, self).addArgs(parser)
        parser.add_argument("--backend", type=str, default="dataframe", help="Backend to use, 'dataframe' (default) or 'lazy'")
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        era = sampleCfg.get("era") if sampleCfg else None
        isMC = self.isMC(sample)
        metName = "METFixEE2017" if era == "2017" else "MET"
        isNotWorker = True # for tests - more realistic: (self.args.distributed != "worker")
        ## Decorate the tree
        from bamboo.treedecorators import NanoAODDescription, nanoRochesterCalc, nanoJetMETCalc, nanoJetMETCalc_METFixEE2017
        tree,noSel,be,lumiArgs = super(NanoZMuMuBase, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg,
                description=NanoAODDescription.get("v5", year=(era if era else "2016"), isMC=isMC,
                    systVariations=[ nanoRochesterCalc, (nanoJetMETCalc_METFixEE2017 if era == "2017" else nanoJetMETCalc) ]), ## will do Jet and MET variations, and the Rochester correction
                lazyBackend=(self.args.backend == "lazy"))
        ## per-year/era options
        puWeightsFile = None
        jecTag, smearTag, jesUncertaintySources = None, None, None
        rochesterFile = None
        if era == "2016":
            rochesterFile = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "data", "RoccoR2016.txt")
            if isMC:
                puWeightsFile = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "data", "puweights.json")
                jecTag = "Summer16_07Aug2017_V20_MC"
                smearTag = "Summer16_25nsV1_MC"
                jesUncertaintySources = ["Total"]
            else:
                if "2016B" in sample or "2016C" in sample or "2016D" in sample:
                    jecTag = "Summer16_07Aug2017BCD_V11_DATA"
                elif "2016E" in sample or "2016F" in sample:
                    jecTag = "Summer16_07Aug2017EF_V11_DATA"
                elif "2016G" in sample or "2016H" in sample:
                    jecTag = "Summer16_07Aug2017GH_V11_DATA"
                else:
                    raise ValueError(f"Could not deduce data JEC tag for sample {sample}")
        elif era == "2017":
            if isMC:
                jecTag = "Fall17_17Nov2017_V32_MC"
                smearTag = "Fall17_V3_MC"
                jesUncertaintySources = ["Total"]
            else:
                if "2017B" in sample:
                    jecTag = "Fall17_17Nov2017B_V32_DATA"
                elif "2017C" in sample:
                    jecTag = "Fall17_17Nov2017C_V32_DATA"
                elif "2017D" in sample or "2017E" in sample:
                    jecTag = "Fall17_17Nov2017DE_V32_DATA"
                elif "2017F" in sample:
                    jecTag = "Fall17_17Nov2017F_V32_DATA"
                else:
                    raise ValueError(f"Could not deduce data JEC tag for sample {sample}")
        ## always-on event weights
        if isMC:
            mcWgts = [ tree.genWeight ]
            if puWeightsFile:
                from bamboo.analysisutils import makePileupWeight
                mcWgts.append(makePileupWeight(puWeightsFile, tree.Pileup_nTrueInt, variation="Nominal",
                    nameHint="bamboo_puWeight{0}".format("".join(c for c in sample if c.isalnum()))))
            else:
                logger.warning("Running on MC without pileup reweighting")
            from bamboo import treefunctions as op
            mcWgts += [
                op.systematic(op.c_float(1.), name="qcdScale", **{f"qcdScale{i:d}" : tree.LHEScaleWeight[i] for i in (0, 1, 3, 5, 7, 8)}),
                op.systematic(op.c_float(1.), name="psISR", up=tree.PSWeight[2], down=tree.PSWeight[0]),
                op.systematic(op.c_float(1.), name="psFSR", up=tree.PSWeight[3], down=tree.PSWeight[1]),
                ]
            noSel = noSel.refine("mcWeight", weight=mcWgts)
        ## configure corrections and variations
        from bamboo.analysisutils import configureJets, configureType1MET, configureRochesterCorrection
        try:
            configureJets(tree._Jet, "AK4PFchs", jec=jecTag, smear=smearTag, jesUncertaintySources=jesUncertaintySources, mayWriteCache=isNotWorker, isMC=isMC, backend=be, uName=sample)
        except Exception as ex:
            logger.exception("Problem while configuring jet correction and variations")
        try:
            configureType1MET(getattr(tree, f"_{metName}"), jec=jecTag, smear=smearTag, jesUncertaintySources=jesUncertaintySources, mayWriteCache=isNotWorker, isMC=isMC, backend=be, uName=sample)
        except Exception as ex:
            logger.exception("Problem while configuring MET correction and variations")
        try:
            configureRochesterCorrection(tree._Muon, rochesterFile, isMC=isMC, backend=be, uName=sample)
        except Exception as ex:
            logger.exception("Problem while configuring the Rochester correction")

        return tree,noSel,be,lumiArgs

class NanoZMuMu(NanoZMuMuBase, NanoAODHistoModule):
    """ Example module: Z->MuMu histograms from NanoAOD """
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import CutFlowReport, SummedPlot, EquidistantBinning
        from bamboo import treefunctions as op
        from bamboo.analysisutils import forceDefine

        era = sampleCfg.get("era") if sampleCfg else None

        plots = []

        ## calculate (corrected) muon 4-momenta before accessing them
        forceDefine(t._Muon.calcProd, noSel)

        muons = op.select(t.Muon, lambda mu : mu.pt > 20.)

        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])
        plots.append(Plot.make1D("dimu_M", op.invariant_mass(muons[0].p4, muons[1].p4), twoMuSel,
                EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

        ## evaluate jet and MET for all events passing twoMuSel
        ## more optimization will be needed with systematics etc.
        forceDefine(t._Jet.calcProd, twoMuSel)
        forceDefine(getattr(t, "_{0}".format("MET" if era != "2017" else "METFixEE2017")).calcProd, twoMuSel)

        jets = op.select(t.Jet, lambda j : j.pt > 20.)
        plots.append(Plot.make1D("nJets", op.rng_len(jets), twoMuSel,
                EquidistantBinning(10, 0., 10.), title="Number of jets"))

        twoMuTwoJetSel = twoMuSel.refine("twoMuonsTwoJets", cut=[ op.rng_len(jets) > 1 ])

        leadjpt = Plot.make1D("leadJetPT", jets[0].pt, twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Leading jet PT")
        subleadjpt = Plot.make1D("subleadJetPT", jets[1].pt, twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="Subleading jet PT")
        plots += [ leadjpt, subleadjpt, SummedPlot("twoLeadJetPT", [leadjpt, subleadjpt], xTitle="Leading two jet PTs") ]
        met = t.MET if era != "2017" else t.METFixEE2017
        plots.append(Plot.make1D("MET", met.pt, twoMuTwoJetSel,
                EquidistantBinning(50, 0., 250.), title="MET PT"))

        selsForCutFlowReport = [ twoMuTwoJetSel ]

        deepCSVFile = os.path.join(os.path.dirname(os.path.dirname(__file__)), "tests", "data", "DeepCSV_2016LegacySF_V1.csv")
        if os.path.exists(deepCSVFile): ## protection for CI tests
            btagReader = makeBTagCalibrationReader("csv", deepCSVFile, wp="Loose", sysType="central", otherSysTypes=("up", "down"), measurementType="comb", sel=noSel, uName=sample)
            #sf_deepcsv = BtagSF(btagReader, getDiscri=lambda j : j.btagDeepB) ## for reshaping
            sf_deepcsv = BtagSF(btagReader)

            bJets_DeepCSVLoose = op.select(jets, lambda j : j.btagDeepB > 0.2217)
            bTagSel = twoMuTwoJetSel.refine("twoMuonsTwoJetsB", cut=[ op.rng_len(bJets_DeepCSVLoose) > 0 ],
                    weight=(sf_deepcsv(bJets_DeepCSVLoose[0]) if self.isMC(sample) else None))
            plots.append(Plot.make1D("bjetpt", bJets_DeepCSVLoose[0].pt, bTagSel, EquidistantBinning(50, 0., 250.), title="B-jet pt"))

            selsForCutFlowReport.append(bTagSel)
        plots.append(CutFlowReport("mumujj", selsForCutFlowReport))

        return plots

class SkimNanoZMuMu(NanoZMuMuBase, NanoAODSkimmerModule):
    def defineSkimSelection(self, tree, noSel, sample=None, sampleCfg=None):
        from bamboo import treefunctions as op
        from bamboo.analysisutils import forceDefine
        forceDefine(tree._Muon.calcProd, noSel)
        muons = op.select(tree.Muon, lambda mu : op.AND(mu.pt > 20., op.abs(mu.eta) < 2.4))
        hasTwoMu = noSel.refine("hasTwoMu", cut=(op.rng_len(muons) >= 2))
        varsToKeep = {"nMuon": None, "Muon_eta": None, "Muon_pt": None} ## from input file
        varsToKeep["nSelMuons"] = op.static_cast("UInt_t", op.rng_len(muons)) ## TBranch doesn't accept size_t
        varsToKeep["selMuons_i"] = muons.idxs
        varsToKeep["selMu_miniPFRelIsoNeu"] = op.map(muons, lambda mu : mu.miniPFRelIso_all - mu.miniPFRelIso_chg)
        return hasTwoMu, varsToKeep
