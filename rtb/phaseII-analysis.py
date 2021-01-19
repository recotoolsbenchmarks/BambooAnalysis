""" Example CMS Phase2 simulation analysis module """

#############################################
##                                         ##
## methods and classes for a common module ##
##                                         ##
#############################################

from bamboo.analysismodules import AnalysisModule, HistogramsModule

class CMSPhase2SimModule(AnalysisModule):
    """ Base module for processing Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimModule, self).__init__(args)
    def prepareTree(self, tree, sample=None, sampleCfg=None):
        from bamboo.treedecorators import decorateCMSPhase2SimTree
        from bamboo.dataframebackend import DataframeBackend
        t = decorateCMSPhase2SimTree(tree, isMC=True)
        be, noSel = DataframeBackend.create(t)
        return t, noSel, be, tuple()
    def mergeCounters(self, outF, infileNames, sample=None):
        from bamboo.root import gbl
        outF.cd()
        hNEvts = gbl.TH1F("nEvents", "Number of events", 1, 0., 1.)
        for fn in infileNames:
            f = gbl.TFile.Open(fn)
            tup = f.Get("myana/mytree") ## FIXME
            if not tup:
                raise RuntimeError("File {0} does not have a tree {1}".format(fn, self.args.treeName))
            hNEvts.Fill(0, tup.GetEntries())
        outF.cd()
        hNEvts.Write("nEvents")
    def readCounters(self, resultsFile):
        hNEvts = resultsFile.Get("nEvents")
        return {"nEvents": hNEvts.GetBinContent(1)}

class CMSPhase2SimHistoModule(CMSPhase2SimModule, HistogramsModule):
    """ Base module for producing plots from Phase2 flat trees """
    def __init__(self, args):
        super(CMSPhase2SimHistoModule, self).__init__(args)


################################
##                            ##
## The actual analysis module ##
##                            ##
################################

class CMSPhase2SimTest(CMSPhase2SimHistoModule):
    """ Example plotter module for Phase2 flat trees """
    def definePlots(self, t, noSel, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, CutFlowReport, EquidistantBinning, VariableBinning
        from bamboo import treefunctions as op

        plots = []

        #definitions

        electrons = op.select(t.elec, lambda el : op.AND(
        el.pt > 20., op.abs(el.eta) < 2.5
        ))
        
        muons = op.select(t.muon, lambda mu : op.AND(
        mu.pt > 20., op.abs(mu.eta) < 2.5
        ))
        
        cleanedElectrons = op.select(electrons, lambda el : op.NOT(
        op.rng_any(muons, lambda mu : op.deltaR(el.p4, mu.p4) < 0.3 )
        ))
        
        isolatedElectrons = op.select(cleanedElectrons, lambda ce : ce.isopass & (1<<2) )
        
        identifiedElectrons = op.select(isolatedElectrons, lambda ie : ie.idpass & (1<<2) )
        
        cleanedMuons = op.select(muons, lambda mu : op.NOT(
        op.rng_any(electrons, lambda el : op.deltaR(mu.p4, el.p4) < 0.3 )
        ))
        
        isolatedMuons = op.select(cleanedMuons, lambda cm : cm.isopass & (1<<2) )
        
        identifiedMuons = op.select(isolatedMuons, lambda im : im.idpass & (1<<2) )
        
        InvMassMuMU = op.invariant_mass(identifiedMuons[0].p4, identifiedMuons[1].p4 )
        
        cleanedJets = op.select(t.jetpuppi, lambda j : op.AND(
        op.NOT(op.rng_any(identifiedElectrons, lambda el : op.deltaR(el.p4, j.p4) < 0.3) ),
        op.NOT(op.rng_any(identifiedMuons, lambda mu : op.deltaR(mu.p4, j.p4) < 0.3) )
        ))
        
        cleanedGoodJets = op.select(cleanedJets, lambda j : op.AND(
        j.pt > 30, op.abs(j.eta) < 2.5
        ))
        
        #selections

        sel1 = noSel.refine("nmumu", cut = [op.AND((op.rng_len(identifiedMuons) > 1), (op.product(identifiedMuons[0].charge, identifiedMuons[1].charge) < 0 ))]) # Oppositely charged MuMu selection

        sel2 = sel1.refine("InvM", cut = [op.NOT(op.in_range(76, InvMassMuMU, 106))]) # Invariant mass selection
            
        sel3 = sel2.refine("njet", cut = [op.rng_len(cleanedGoodJets) > 1]) # two jets selection
            
        sel4 = sel3.refine("btag", cut = [op.OR(cleanedGoodJets[0].btag & (1<<2), cleanedGoodJets[1].btag & (1<<2))]) # at least 1 among two leading jets is b-tagged, we are taking the second isopass to be on which is equal to the medium b-tag working point
            
        sel5 = sel4.refine("MET", cut = [t.metpf[0].pt > 40]) # MET > 40 GeV, here the met pt cut requires indexing and t.met[0].pt is equivalent to t.met.pt
                    
        #plots
            
            #noSel
        plots.append(Plot.make1D("nJetsNoSel", op.rng_len(cleanedGoodJets), noSel, EquidistantBinning(10, 0., 10.), title="nJets"))
            
        plots.append(Plot.make1D("nMuNoSel", op.rng_len(identifiedMuons), noSel, EquidistantBinning(15, 0., 15.), title="nMuons"))
        
        plots.append(Plot.make1D("METptNoSel", t.metpf[0].pt, noSel, EquidistantBinning(50, 0., 250), title="MET_PT"))
        
            #sel1
        
        plots.append(Plot.make1D("nJetsSel1", op.rng_len(cleanedGoodJets), sel1, EquidistantBinning(10, 0., 10.), title="nJets"))  
            
        plots.append(Plot.make1D("nMuSel1", op.rng_len(identifiedMuons), sel1, EquidistantBinning(10, 0., 10.), title="nMuMUbar PT > 20, |eta| < 2.5"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel1", InvMassMuMU, sel1, EquidistantBinning(30, 0, 300), title="m(ll) < 76 & m(ll) > 106 GeV"))
                  
        plots.append(Plot.make1D("LeadingMuonPTSel1", muons[0].pt, sel1, EquidistantBinning(30, 0., 250.), title=" Leading Muon PT"))
            
        plots.append(Plot.make1D("SubLeadingMuonPTSel1", muons[1].pt, sel1, EquidistantBinning(30, 0., 250.), title="SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel1", muons[0].eta, sel1, EquidistantBinning(30, -3, 3), title=" Leading Muon eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel1", muons[1].eta, sel1, EquidistantBinning(30, -3, 3), title="SubLeading Muon eta"))
        
        plots.append(Plot.make1D("METptSel1", t.metpf[0].pt, sel1, EquidistantBinning(50, 0., 250), title="MET_PT"))
                            
            #sel2

        plots.append(Plot.make1D("nJetsSel2", op.rng_len(cleanedGoodJets), sel2, EquidistantBinning(10, 0., 10.), title="nJets"))

        plots.append(Plot.make1D("nMuSel2", op.rng_len(identifiedMuons), sel2, EquidistantBinning(10, 0., 10.), title="nMuMUbar PT > 20, |eta| < 2.5"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel2", InvMassMuMU, sel2, EquidistantBinning(20, 20., 300.), title="m(ll) < 76 & m(ll) > 106 GeV"))
        
        plots.append(Plot.make1D("LeadingMuonPTSel2", muons[0].pt, sel2, EquidistantBinning(30, 0., 250.), title=" Leading Muon PT"))

        plots.append(Plot.make1D("SubLeadingMuonPTSel2", muons[1].pt, sel2, EquidistantBinning(30, 0., 200.), title=" SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel2", muons[0].eta, sel2, EquidistantBinning(30, -3, 3), title=" Leading Muon Eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel2", muons[1].eta, sel2, EquidistantBinning(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("METptSel2", t.metpf[0].pt, sel2, EquidistantBinning(50, 0., 250), title="MET_PT"))
                        
            #sel3
            
        plots.append(Plot.make1D("nJetsSel3", op.rng_len(cleanedGoodJets), sel3, EquidistantBinning(10, 0., 10.), title="nJets > 2"))

        plots.append(Plot.make1D("LeadingJetPTSel3", cleanedGoodJets[0].pt, sel3, EquidistantBinning(50, 0., 350.), title="Leading jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel3", cleanedGoodJets[1].pt, sel3, EquidistantBinning(50, 0., 350.), title="SubLeading jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel3", cleanedGoodJets[0].eta, sel3, EquidistantBinning(30, -3, 3), title="Leading jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel3", cleanedGoodJets[1].eta, sel3, EquidistantBinning(30, -3, 3), title="SubLeading jet Eta"))

        plots.append(Plot.make1D("nMuSel3", op.rng_len(identifiedMuons), sel3, EquidistantBinning(10, 0., 10.), title="nMuMUbar PT > 20, |eta| < 2.5"))    
        
        plots.append(Plot.make1D("LeadingMuonPTSel3", muons[0].pt, sel3, EquidistantBinning(30, 0., 250.), title=" Leading Muon PT"))
            
        plots.append(Plot.make1D("SubLeadingMuonPTSel3", muons[1].pt, sel3, EquidistantBinning(30, 0., 200.), title=" SubLeading Muon PT"))

        plots.append(Plot.make1D("LeadingMuonEtaSel3", muons[0].eta, sel3, EquidistantBinning(30, -3, 3), title=" Leading Muon Eta"))
            
        plots.append(Plot.make1D("SubLeadingMuonEtaSel3", muons[1].eta, sel3, EquidistantBinning(30, -3, 3), title=" SubLeading Muon Eta"))
                
        plots.append(Plot.make1D("InvMassTwoMuonsSel3", InvMassMuMU, sel3, EquidistantBinning(30, 0, 300), title="m(ll) < 76 & m(ll) > 106 GeV"))
        
        plots.append(Plot.make1D("METptSel3", t.metpf[0].pt, sel3, EquidistantBinning(50, 0., 250), title="MET_PT"))
        
            #sel4
             
        plots.append(Plot.make1D("nJetsSel4", op.rng_len(cleanedGoodJets), sel4, EquidistantBinning(10, 0, 10), title="B-tagged jets"))
            
        plots.append(Plot.make1D("LeadingJetPTSel4", cleanedGoodJets[0].pt, sel4, EquidistantBinning(50, 0., 250.), title="Leading b-tagged jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel4", cleanedGoodJets[1].pt, sel4, EquidistantBinning(50, 0., 250.), title="SubLeading b-tagged jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel4", cleanedGoodJets[0].eta, sel4, EquidistantBinning(30, -3, 3.), title="Leading b-tagged jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel4", cleanedGoodJets[1].eta, sel4, EquidistantBinning(30, -3, 3.), title="SubLeading b-tagged jet Eta"))
        
        plots.append(Plot.make1D("nMuSel4", op.rng_len(identifiedMuons), sel4, EquidistantBinning(10, 0., 10.), title="nMuMUbar PT > 20, |eta| < 2.5"))
                
        plots.append(Plot.make1D("LeadingMuonPTSel4", muons[0].pt, sel4, EquidistantBinning(30, 0., 250.), title=" Leading Muon PT"))
        
        plots.append(Plot.make1D("SubLeadingMuonPTSel4", muons[1].pt, sel4, EquidistantBinning(30, 0., 200.), title=" SubLeading Muon PT"))
        
        plots.append(Plot.make1D("LeadingMuonEtaSel4", muons[0].eta, sel4, EquidistantBinning(30, -3, 3), title=" Leading Muon Eta"))
        
        plots.append(Plot.make1D("SubLeadingMuonEtaSel4", muons[1].eta, sel4, EquidistantBinning(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel4", InvMassMuMU, sel4, EquidistantBinning(30, 0, 300), title="m(ll) < 76 & m(ll) > 106 GeV"))
        
        plots.append(Plot.make1D("METptSel4", t.metpf[0].pt, sel4, EquidistantBinning(50, 0., 250), title="MET_PT"))

            #sel5
                
        plots.append(Plot.make1D("nJetsSel5", op.rng_len(cleanedGoodJets), sel5, EquidistantBinning(10, 0, 10), title="B-tagged jets"))
            
        plots.append(Plot.make1D("LeadingJetPTSel5", cleanedGoodJets[0].pt, sel5, EquidistantBinning(50, 0., 250.), title="Leading b-tagged jet PT"))
            
        plots.append(Plot.make1D("SubLeadingJetPTSel5", cleanedGoodJets[1].pt, sel5, EquidistantBinning(50, 0., 250.), title="SubLeading b-tagged jet PT"))
            
        plots.append(Plot.make1D("LeadingJetEtaSel5", cleanedGoodJets[0].eta, sel5, EquidistantBinning(30, -3, 3.), title="Leading b-tagged jet Eta"))
            
        plots.append(Plot.make1D("SubLeadingJetEtaSel5", cleanedGoodJets[1].eta, sel5, EquidistantBinning(30, -3, 3.), title="SubLeading b-tagged jet Eta"))
        
        plots.append(Plot.make1D("nMuSel5", op.rng_len(identifiedMuons), sel5, EquidistantBinning(10, 0., 10.), title="nMuMUbar PT > 20, |eta| < 2.5"))
                
        plots.append(Plot.make1D("LeadingMuonPTSel5", muons[0].pt, sel5, EquidistantBinning(30, 0., 250.), title=" Leading Muon PT"))
        
        plots.append(Plot.make1D("SubLeadingMuonPTSel5", muons[1].pt, sel5, EquidistantBinning(30, 0., 200.), title=" SubLeading Muon PT"))
        
        plots.append(Plot.make1D("LeadingMuonEtaSel5", muons[0].eta, sel5, EquidistantBinning(30, -3, 3), title=" Leading Muon Eta"))
        
        plots.append(Plot.make1D("SubLeadingMuonEtaSel5", muons[1].eta, sel5, EquidistantBinning(30, -3, 3), title=" SubLeading Muon Eta"))

        plots.append(Plot.make1D("InvMassTwoMuonsSel5", InvMassMuMU, sel5, EquidistantBinning(30, 0, 300), title="m(ll) < 76 & m(ll) > 106 GeV"))

        plots.append(Plot.make1D("METptSel5", t.metpf[0].pt, sel5, EquidistantBinning(50, 0., 250), title="MET_PT > 40"))

        # Efficiency Report on terminal and the .tex output

        cfr = CutFlowReport("yields")
        cfr.add(noSel, "No selection")
        cfr.add(sel1, "nMuMu >= 2")
        cfr.add(sel2, "InvM")
        cfr.add(sel3, "nJet >= 2")
        cfr.add(sel4, "btag")
        cfr.add(sel5, "MET")

        plots.append(cfr)
                            
        return plots
