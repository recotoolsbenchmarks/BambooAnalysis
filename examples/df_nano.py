# A first attempt at `ROOT::RDataFrame` wrappers in python
# ========================================================
#
# *Warning*: this is very much work in progress... the idea is to wrap a NanoAOD in an analyst-friendly way,
# to make writing analysis code (from tree to plots) as flexible as possible,
# such that python scripting can be used to avoid repetition and keep the code compact and maintainable.
#
# A few imports to get started, and debug logging such that we can see what is going on.
import logging
logging.basicConfig(level=logging.DEBUG)
from bamboo.treedecorators import decorateNanoAOD
from bamboo.dataframebackend import DataframeBackend
from bamboo.analysisutils import addLumiMask
import ROOT
if hasattr(ROOT, "enableJSVis"): ROOT.enableJSVis() ## for inline plots in the notebook
# `decorateNanoAOD` will wrap the NanoAOD in a more pythonic set of objects, `DataframeBackend` takes care of
# passing the right instructions to a `ROOT::RDataFrame` instance when we define plots.
f = ROOT.TFile.Open("NanoAOD_SingleMu_test.root")
t = decorateNanoAOD(f.Get("Events"))
be, noSel = DataframeBackend.create(t)
noSel = addLumiMask(noSel, "Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt",
            runRange=(276831, 277420), runAndLS=(t.run, t.luminosityBlock))
# At this point we have a wrapped tree `t`, a reference to the backend
# (which we will only need in the end, to retrieve the plot results), 
# and a basic selection `noSel`, which lets all events pass.
# To define a plot, we will need to pass (at the very least) a name for it,
# a variable to plot (with a suitable binning), and a selection to apply
# (which could be noSel if the plot is to be filled for all records).
from bamboo.plots import Plot, EquidistantBinning
import bamboo.treefunctions as op
# The `Plot` class is defined in the plots module, together with helper classes
# (e.g. binnings, such as `EquidistantBinning` and `VariableBinning`).
# The `treefunctions` module defines mathematical functions etc. for the
# tree proxy objects.

# Let's start with an example to demonsrate the power of the plot wrappers:
# a dimuon mass plot. First we'll need a dimuon selection
twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(t.Muon) > 1 ])
# Then we can define the plot, using the momenta of the two leading muons.
# Note that the `p4` branch (or `Muon_p4`) is not in the original tree.
# It is transparently created on-the-fly from the pt, eta, phi and mass
# branches (the backend will tell the `RDataFrame` to construct it when needed).
myMPlot = Plot.make1D("dimu_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSel,
            EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass")
# The original pt attribute of the muon (from `Muon_pt`) is als still there:
myPlot = Plot.make1D("dimu_Pt2", t.Muon[1].pt, twoMuSel,
            EquidistantBinning(50,0., 200.), title="Trailing muon PT")
myOtherPlot = Plot.make1D("dimu_Pt1", t.Muon[0].p4.Pt(), twoMuSel,
            EquidistantBinning(50, 0., 200.), title="Leading muon PT")

# After this warmup, let's get a bit more realistic, and define a selected
# muon collection and plot the pT of the first one (if there are).
myMu = op.select(t.Muon, lambda mu : op.abs(mu.eta) < 2.)
myNsPlot = Plot.make1D("nMu10", op.rng_len(myMu), noSel,
            EquidistantBinning(10, 0., 10.), title="Number of muons with p_{T} > 10 GeV")
hasMyMuSel = noSel.refine("oneMyMuon", cut=[ op.rng_len(myMu) > 0 ])
leadMyMuPtPlot = Plot.make1D("leadMymuPt", myMu[0].pt, hasMyMuSel,
            EquidistantBinning(50, 0., 200.), title="Leading myMuon PT")
# If we only needed one central muon, we could also have used `rng_any` and `rng_find`:
isVeryCentral = lambda mu : op.abs(mu.eta) < 2.
hasCentralMuon = noSel.refine("hasCentralMu", cut=[op.rng_any(t.Muon, isVeryCentral)])
leadCentralMu = op.rng_find(t.Muon, isVeryCentral)
centralMuPtPlot = Plot.make1D("centMu_Pt1", leadCentralMu.pt, hasCentralMuon,
            EquidistantBinning(50, 0., 200.), title="Leading central Muon PT")
# Range functions that (unlike `rng_find` and `rng_any`) always process
# the complete range are also defined, let's for example sum the pseudorapidities
# of the very central, and then all muons (despite the limited physical interest):
myMuSumEtaPlot = Plot.make1D("sumMyMuEta", op.rng_sum(myMu, lambda mu : mu.eta), noSel,
            EquidistantBinning(50, -10., 10.), title="Sum of myMu eta")
muSumEtaPlot = Plot.make1D("sumMuEta", op.rng_sum(t.Muon, lambda mu : mu.eta), noSel,
            EquidistantBinning(50, -10., 10.), title="Sum of mu eta")
# A very useful function looks for the object with the largest/smallest
# value for some function, `rng_max`/`min_element`, e.g. the most central muon:
mostCentralMuon = op.rng_min_element_by(t.Muon, lambda mu : op.abs(mu.eta))
mostCentralMuPt = Plot.make1D("mostCentralMuPt", mostCentralMuon.pt, hasCentralMuon,
            EquidistantBinning(50, 0., 200.), title="Most central muon PT")

# Also scalefactors are easy to add, as weight factors for the Selection:
from bamboo.scalefactors import get_scalefactor
lSF = get_scalefactor("lepton", ("muon_2016_80", "iso_loose_id_loose"), combine="weight")
twoMuSelSF = noSel.refine("twoMuonsSF", cut=[ op.rng_len(t.Muon) > 1 ], weight=[ lSF(t.Muon[0], withMCCheck=False), lSF(t.Muon[1], withMCCheck=False) ])
newMPlot = Plot.make1D("dimuSF_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSelSF,
            EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass")

# The basic infrastructure for on-the-fly jet systematic variations is also in place,
# this example does not modify the jet branches, but could be made to do so by calling
# a method in C++. `t.Jet.nominal` (or `t.Jet["nominal"]`) has the original jet,
# `t.Jet[variation]` a modified jet container.
op.addKinematicVariation(t.Jet, "no-op", modif=( lambda j : j.p4 ), pred=( lambda mp4, j : mp4.Pt() > 10. ))
nJets = Plot.make1D("nJets", op.rng_len(t.Jet["no-op"]), noSel, EquidistantBinning(10, 0., 10.), title="Unmodified jet multiplicity")

# Making two-particle combinations (higher N to be added) is also easy:
osdimu = op.combine(t.Muon, N=2, pred=lambda m1, m2 : m1.charge != m2.charge)
oselmu = op.combine((t.Electron, t.Muon), pred=lambda el, mu : el.charge != mu.charge)
nOSDimu = Plot.make1D("nOsDimu", op.rng_len(osdimu), noSel, EquidistantBinning(10, 0., 10.), title="Number of opposite-sign dimuon pairs")
nOSelmu = Plot.make1D("nOselmu", op.rng_len(oselmu), noSel, EquidistantBinning(10, 0., 10.), title="Number of opposite-sign electron-muon pairs")
## TODO stress-test: use selected muons, select afterwards, combine electrons and muons

tagEle   = op.select(t.Electron, lambda el : op.AND(el.pt >= 30., op.abs(el.eta+el.deltaEtaSC) < 2.1, op.NOT(op.in_range(1.4442, op.abs(el.eta+el.deltaEtaSC), 1.566)),
    op.rng_any(t.TrigObj, lambda to : op.AND(to.id == 11, to.filterBits & 2, op.sqrt((to.eta-el.eta)**2+(to.phi-el.phi)**2) < .4)))) ## TODO review 'filterBits'
nTagEle = Plot.make1D("nTagEle", op.rng_len(tagEle), noSel, EquidistantBinning(10, 0., 10.), title="Number of tag electrons (partial selection)")
probeEle = op.select(t.Electron, lambda el : op.abs(el.eta+el.deltaEtaSC) < 2.5)
eleTagAndProbe = op.combine((tagEle, probeEle), pred=lambda elt, elp : op.in_range(50., op.invariant_mass(elt.p4, elp.p4), 130.), samePred=lambda o1,o2: o1.idx != o2.idx)
hasTagAndProbeEle = noSel.refine("hasTagAndProbe", cut=(op.rng_len(eleTagAndProbe) > 0)) ## TODO add trigger event selection
theTnPpair = eleTagAndProbe[0] ## TODO add more strategies
tnpelemass = Plot.make1D("eleTnPmassInteg", op.invariant_mass(theTnPpair[0].p4, theTnPpair[1].p4), hasTagAndProbeEle, EquidistantBinning(80, 50., 130.), title="Tag&Probe electron pair invariant mass")

# That's a nice start - but we haven't plotted any histograms yet!
# The `ROOT::RDataFrame` will do that as soon as we ask for one of them:
cv = ROOT.TCanvas("c1", "A few dimuon plots", 1500, 350)
cv.Divide(3)
cv.cd(1)
be.getResults(myPlot)[0].Draw()
cv.cd(2)
be.getResults(myOtherPlot)[0].Draw()
cv.cd(3)
be.getResults(myMPlot)[0].Draw()
cv.Update()
cv.Draw()

cv2 = ROOT.TCanvas("c2")
cv2.Divide(3)
cv2.cd(1)
be.getResults(myNsPlot)[0].Draw()
cv2.cd(2)
be.getResults(leadMyMuPtPlot)[0].Draw()
cv2.cd(3)
be.getResults(centralMuPtPlot)[0].Draw()
cv2.Update()
cv2.Draw()

cv3 = ROOT.TCanvas("c3")
cv3.Divide(2)
cv3.cd(1)
be.getResults(myMuSumEtaPlot)[0].Draw()
cv3.cd(2)
be.getResults(mostCentralMuPt)[0].Draw()
cv3.Update()
cv3.Draw()

cv4 = ROOT.TCanvas("c4")
cv4.Divide(2)
cv4.cd(1)
be.getResults(newMPlot)[0].Draw()
cv4.cd(2)
be.getResults(nJets)[0].Draw()
cv4.Update()

cv5 = ROOT.TCanvas("c5")
cv5.Divide(2)
cv5.cd(1)
be.getResults(nOSDimu)[0].Draw()
cv5.cd(2)
be.getResults(nOSelmu)[0].Draw()
cv5.Update()

cv6 = ROOT.TCanvas("c6")
cv6.Divide(2)
cv6.cd(1)
be.getResults(nTagEle)[0].Draw()
cv6.cd(1)
be.getResults(tnpelemass)[0].Draw()
cv6.Update()

# There are a few more things to figure out:
#   - integrate with slurm. The easiest is probably to have a common module
#     that runs either in a single thread or on slurm, that gets a list of plots
#     from a method `definePlots(tree, noSel)`, and saves all histograms to a file.
#   - Combinatorics (pairs and triplets of objects)
#   - "skims" (should be supported through the `Snapshot` method)
#   - other tree formats?
#   - `your idea here`
