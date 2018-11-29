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
import ROOT
if hasattr(ROOT, "enableJSVis"): ROOT.enableJSVis() ## for inline plots in the notebook
# `decorateNanoAOD` will wrap the NanoAOD in a more pythonic set of objects, `DataframeBackend` takes care of
# passing the right instructions to a `ROOT::RDataFrame` instance when we define plots.
f = ROOT.TFile.Open("examples/NanoAOD_SingleMu_test.root")
t = decorateNanoAOD(f.Get("Events"))
be, noSel = DataframeBackend.create(t)
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

# That's a nice start - but we haven't plotted any histograms yet!
# The `ROOT::RDataFrame` will do that as soon as we ask for one of them:
cv = ROOT.TCanvas("c1", "A few dimuon plots", 1500, 350)
cv.Divide(3)
cv.cd(1)
be.getPlotResult(myPlot).Draw()
cv.cd(2)
be.getPlotResult(myOtherPlot).Draw()
cv.cd(3)
be.getPlotResult(myMPlot).Draw()
cv.Update()
cv.Draw()

cv2 = ROOT.TCanvas("c2")
cv2.Divide(3)
cv2.cd(1)
be.getPlotResult(myNsPlot).Draw()
cv2.cd(2)
be.getPlotResult(leadMyMuPtPlot).Draw()
cv2.cd(3)
be.getPlotResult(centralMuPtPlot).Draw()
cv2.Update()
cv2.Draw()

cv3 = ROOT.TCanvas("c3")
cv3.Divide(2)
cv3.cd(1)
be.getPlotResult(myMuSumEtaPlot).Draw()
cv3.cd(2)
be.getPlotResult(mostCentralMuPt).Draw()
cv3.Update()
cv3.Draw()

# There are a few more things to figure out:
#   - integrate with slurm. The easiest is probably to have a common module
#     that runs either in a single thread or on slurm, that gets a list of plots
#     from a method `definePlots(tree, noSel)`, and saves all histograms to a file.
#   - scale factors: essentially solved, need to copy the code
#   - systematics (I have an idea for jet systematics, need to implement),
#     the interface could become `definePlots(tree, noSel, systVar=None)`, such that
#     the user can choose the right scalefactor or jet transformation based on `systVar`,
#     but without changing the rest of the code.
#   - Combinatorics (pairs and triplets of objects)
#   - "skims" (should be supported through the `Snapshot` method)
#   - other tree formats?
#   - a cool package name
#   - `your idea here`
