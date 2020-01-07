"""
Example analysis module: electron Tag&Probe from a NanoAOD
"""
from bamboo.analysismodules import NanoAODHistoModule
import logging
logger = logging.getLogger(__name__)

class TagAndProbe1(NanoAODHistoModule):
    def __init__(self, args):
        super(TagAndProbe1, self).__init__(args)

    def make1D(self, name, variable, selection, binning, **kwargs):
        from bamboo.plots import Plot
        nm = "{0}_{1}".format(name, selection.name)
        return Plot.make1D(nm, variable, selection, binning, **kwargs)

    def definePlots(self, t, noSel, era=None, sample=None, sampleCfg=None):
        from bamboo.plots import Plot, EquidistantBinning
        from bamboo import treefunctions as op
        from itertools import count

        plots = []

        tnpTrigSel = noSel.refine("hasTnPTrig", cut=t.HLT.Ele32_eta2p1_WPTight_Gsf)

        tagEle = op.select(t.Electron, lambda el : op.AND(el.pt >= 30.,
            op.abs(el.eta+el.deltaEtaSC) < 2.1, op.NOT(op.in_range(1.4442, op.abs(el.eta+el.deltaEtaSC), 1.566)),
            el.cutBased_Fall17_V1,
            op.rng_any(t.TrigObj, lambda to : op.AND(
                to.id == 11, to.filterBits & 2, ## TODO review 'filterBits'
                op.sqrt((to.eta-el.eta)**2+(to.phi-el.phi)**2) < .4
                ))
            ))

        probeEle = op.select(t.Electron, lambda el : op.abs(el.eta+el.deltaEtaSC) < 2.5)

        plots.append(Plot.make1D("nTagEle", op.rng_len(tagEle), tnpTrigSel,
            EquidistantBinning(10, 0., 10.), title="Number of tag electrons (partial selection)"))
        plots.append(Plot.make1D("nProbeEle", op.rng_len(probeEle), tnpTrigSel,
            EquidistantBinning(10, 0., 10.), title="Number of probe electrons (partial selection)"))

        eleTagAndProbe = op.combine((tagEle, probeEle),
                pred=lambda elt, elp : op.in_range(50., op.invariant_mass(elt.p4, elp.p4), 130.),
                samePred=lambda o1,o2: o1.idx != o2.idx)

        plots.append(Plot.make1D("nTnPElePairs", op.rng_len(eleTagAndProbe), tnpTrigSel,
            EquidistantBinning(10, 0., 10.), title="Tag&Probe electron: number of pairs"))

        hasTnPpair = tnpTrigSel.refine("hasTnPpair", cut=(op.rng_len(eleTagAndProbe) > 0))
        theTnPpair = eleTagAndProbe[0] ## TODO add more strategies
        ## theTnpPair = op.rng_pickRandom(tnpPairs, seed=...)
        tag, probe = theTnPpair[0], theTnPpair[1]
        tnpMll = op.invariant_mass(theTnPpair[0].p4, theTnPpair[1].p4)
        plots.append(Plot.make1D("eleTnP_ptTag", theTnPpair[0].p4.Pt(), hasTnPpair,
            EquidistantBinning(50, 0., 150.), title="Tag&Probe electron pair: probe PT"))
        plots.append(Plot.make1D("eleTnP_ptProbe", theTnPpair[1].p4.Pt(), hasTnPpair,
            EquidistantBinning(50, 0., 150.), title="Tag&Probe electron pair: probe PT"))

        plots.append(Plot.make1D("eleTnPmassInteg", tnpMll, hasTnPpair,
            EquidistantBinning(80, 50., 130.), title="Tag&Probe electron pair invariant mass"))

        ptBinning = [ 15., 20., 25., 30., 40., 50., 70., 100., 200. ]
        denomSels = dict(("PT{0}".format(i), hasTnPpair.refine("eleTnP_denom_PT{0}".format(i),
                cut=op.in_range(ptMn, probe.p4.Pt(), ptMx)))
            for i, ptMn, ptMx in zip(count(), ptBinning[:-1], ptBinning[1:]))
        for dSel in denomSels.values():
            plots.append(Plot.make1D("{0}_mll".format(dSel.name), tnpMll, dSel,
                EquidistantBinning(80, 50., 130.), title="Tag&Probe electron pair invariant mass"))

        selsToTest = {
              #"MVAIso90" : op.construct("bool", (probe.mvaFall17V1Iso_WP90,)) ## reading vector<bool> branches needs a fix from ROOT 6.16; 'op.construct' forces copy to bool, this workaround should not be needed (if the bool comes from an expression it isn't)
              "mvaTTH090" : probe.mvaTTH > 0.95 ## for demonstration only
              }
        for tsName, testSel in selsToTest.items():
            for bNm,dSel in denomSels.items():
                nomSel = dSel.refine("eleTnP_pass{0}_{1}".format(tsName, bNm), cut=testSel)
                plots.append(Plot.make1D("{0}_mll".format(nomSel.name), tnpMll, nomSel,
                    EquidistantBinning(80, 50., 130.), title="Tag&Probe electron pair invariant mass"))

        return plots

    ## TODO in postprocess: fit for the Z->e+e- yields and make ratios
