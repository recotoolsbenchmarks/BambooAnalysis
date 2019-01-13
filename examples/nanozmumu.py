"""
Example analysis module: make a dimuon mass plot from a NanoAOD
"""
from bamboo.analysismodules import NanoAODHistoModule

class NanoZMuMu(NanoAODHistoModule):
    """ Example module: Z->MuMu histograms from NanoAOD """
    def __init__(self, args):
        super(NanoZMuMu, self).__init__(args)

    def definePlots(self, t, noSel, systVar="nominal"):
        from bamboo.plots import Plot, EquidistantBinning
        from bamboo import treefunctions as op

        plots = []

        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(t.Muon) > 1 ])
        plots.append(Plot.make1D("dimu_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSel,
                EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

        return plots
