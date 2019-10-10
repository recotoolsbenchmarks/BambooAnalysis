#!/usr/bin/env python
"""
A script to generate a BinnedValues-JSON file for pileup reweighting of MC
"""
import json
import logging
logger = logging.getLogger(__name__)
import numpy as np

mcPUProfiles = {
        #========#
        #  2018  #
        #========#
        ## mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2018_25ns_JuneProjectionFull18_PoissonOOTPU_cfi.py#L3-L32
        "Autumn18_25ns": (
            np.linspace(0., 100., 101),
            [4.695341e-10, 1.206213e-06, 1.162593e-06, 6.118058e-06, 1.626767e-05, 3.508135e-05, 7.12608e-05, 0.0001400641, 0.0002663403, 0.0004867473, 0.0008469, 0.001394142, 0.002169081, 0.003198514, 0.004491138, 0.006036423, 0.007806509, 0.00976048, 0.0118498, 0.01402411, 0.01623639, 0.01844593, 0.02061956, 0.02273221, 0.02476554, 0.02670494, 0.02853662, 0.03024538, 0.03181323, 0.03321895, 0.03443884, 0.035448, 0.03622242, 0.03674106, 0.0369877, 0.03695224, 0.03663157, 0.03602986, 0.03515857, 0.03403612, 0.0326868, 0.03113936, 0.02942582, 0.02757999, 0.02563551, 0.02362497, 0.02158003, 0.01953143, 0.01750863, 0.01553934, 0.01364905, 0.01186035, 0.01019246, 0.008660705, 0.007275915, 0.006043917, 0.004965276, 0.004035611, 0.003246373, 0.002585932, 0.002040746, 0.001596402, 0.001238498, 0.0009533139, 0.0007282885, 0.000552306, 0.0004158005, 0.0003107302, 0.0002304612, 0.0001696012, 0.0001238161, 8.96531e-05, 6.438087e-05, 4.585302e-05, 3.23949e-05, 2.271048e-05, 1.580622e-05, 1.09286e-05, 7.512748e-06, 5.140304e-06, 3.505254e-06, 2.386437e-06, 1.625859e-06, 1.111865e-06, 7.663272e-07, 5.350694e-07, 3.808318e-07, 2.781785e-07, 2.098661e-07, 1.642811e-07, 1.312835e-07, 1.081326e-07, 9.141993e-08, 7.890983e-08, 6.91468e-08, 6.119019e-08, 5.443693e-08, 4.85036e-08, 4.31486e-08, 3.822112e-08]
            ),
        #========#
        #  2017  #
        #========#
        ## mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2017_25ns_WinterMC_PUScenarioV1_PoissonOOTPU_cfi.py#L3-L112
        "Fall17_25ns": (
            np.linspace(0., 99., 100),
            [3.39597497605e-05, 6.63688402133e-06, 1.39533611284e-05, 3.64963078209e-05, 6.00872171664e-05, 9.33932578027e-05, 0.000120591524486, 0.000128694546198, 0.000361697233219, 0.000361796847553, 0.000702474896113, 0.00133766053707, 0.00237817050805, 0.00389825605651, 0.00594546732588, 0.00856825906255, 0.0116627396044, 0.0148793350787, 0.0179897368379, 0.0208723871946, 0.0232564170641, 0.0249826433945, 0.0262245860346, 0.0272704617569, 0.0283301107549, 0.0294006137386, 0.0303026836965, 0.0309692426278, 0.0308818046328, 0.0310566806228, 0.0309692426278, 0.0310566806228, 0.0310566806228, 0.0310566806228, 0.0307696426944, 0.0300103336052, 0.0288355370103, 0.0273233309106, 0.0264343533951, 0.0255453758796, 0.0235877272306, 0.0215627588047, 0.0195825559393, 0.0177296309658, 0.0160560731931, 0.0146022004183, 0.0134080690078, 0.0129586991411, 0.0125093292745, 0.0124360740539, 0.0123547104433, 0.0123953922486, 0.0124360740539, 0.0124360740539, 0.0123547104433, 0.0124360740539, 0.0123387597772, 0.0122414455005, 0.011705203844, 0.0108187105305, 0.00963985508986, 0.00827210065136, 0.00683770076341, 0.00545237697118, 0.00420456901556, 0.00367513566191, 0.00314570230825, 0.0022917978982, 0.00163221454973, 0.00114065309494, 0.000784838366118, 0.000533204105387, 0.000358474034915, 0.000238881117601, 0.0001984254989, 0.000157969880198, 0.00010375646169, 6.77366175538e-05, 4.39850477645e-05, 2.84298066026e-05, 1.83041729561e-05, 1.17473542058e-05, 7.51982735129e-06, 6.16160108867e-06, 4.80337482605e-06, 3.06235473369e-06, 1.94863396999e-06, 1.23726800704e-06, 7.83538083774e-07, 4.94602064224e-07, 3.10989480331e-07, 1.94628487765e-07, 1.57888581037e-07, 1.2114867431e-07, 7.49518929908e-08, 4.6060444984e-08, 2.81008884326e-08, 1.70121486128e-08, 1.02159894812e-08]
            ),
        #========#
        #  2016  #
        #========#
        ## mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_Moriond17MC_PoissonOOTPU_cfi.py#L24-L25
        "Moriond17_25ns": (
            np.linspace(0., 75., 76),
	    [1.78653e-05 ,2.56602e-05 ,5.27857e-05 ,8.88954e-05 ,0.000109362 ,0.000140973 ,0.000240998 ,0.00071209 ,0.00130121 ,0.00245255 ,0.00502589 ,0.00919534 ,0.0146697 ,0.0204126 ,0.0267586 ,0.0337697 ,0.0401478 ,0.0450159 ,0.0490577 ,0.0524855 ,0.0548159 ,0.0559937 ,0.0554468 ,0.0537687 ,0.0512055 ,0.0476713 ,0.0435312 ,0.0393107 ,0.0349812 ,0.0307413 ,0.0272425 ,0.0237115 ,0.0208329 ,0.0182459 ,0.0160712 ,0.0142498 ,0.012804 ,0.011571 ,0.010547 ,0.00959489 ,0.00891718 ,0.00829292 ,0.0076195 ,0.0069806 ,0.0062025 ,0.00546581 ,0.00484127 ,0.00407168 ,0.00337681 ,0.00269893 ,0.00212473 ,0.00160208 ,0.00117884 ,0.000859662 ,0.000569085 ,0.000365431 ,0.000243565 ,0.00015688 ,9.88128e-05 ,6.53783e-05 ,3.73924e-05 ,2.61382e-05 ,2.0307e-05 ,1.73032e-05 ,1.435e-05 ,1.36486e-05 ,1.35555e-05 ,1.37491e-05 ,1.34255e-05 ,1.33987e-05 ,1.34061e-05 ,1.34211e-05 ,1.34177e-05 ,1.32959e-05 ,1.33287e-05]
            ),
        ## mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2016_25ns_SpringMC_PUScenarioV1_PoissonOOTPU_cfi.py#L24-L75
        "Spring16_25ns": (
            np.linspace(0., 50., 51),
            [ 0.000829312873542, 0.00124276120498, 0.00339329181587, 0.00408224735376, 0.00383036590008, 0.00659159288946, 0.00816022734493, 0.00943640833116, 0.0137777376066, 0.017059392038, 0.0213193035468, 0.0247343174676, 0.0280848773878, 0.0323308476564, 0.0370394341409, 0.0456917721191, 0.0558762890594, 0.0576956187107, 0.0625325287017, 0.0591603758776, 0.0656650815128, 0.0678329011676, 0.0625142146389, 0.0548068448797, 0.0503893295063, 0.040209818868, 0.0374446988111, 0.0299661572042, 0.0272024759921, 0.0219328403791, 0.0179586571619, 0.0142926728247, 0.00839941654725, 0.00522366397213, 0.00224457976761, 0.000779274977993, 0.000197066585944, 7.16031761328e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
            ),
        #========#
        #  2015  #
        #========#
        ## mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2015_25ns_FallMC_matchData_PoissonOOTPU_cfi.py#L24-L75
        "Fall15_25ns": (
            np.linspace(0., 50., 51),
            [0.000108643, 0.000388957, 0.000332882, 0.00038397, 0.000549167, 0.00105412, 0.00459007, 0.0210314, 0.0573688, 0.103986, 0.142369, 0.157729, 0.147685, 0.121027, 0.08855, 0.0582866, 0.0348526, 0.019457, 0.0107907, 0.00654313, 0.00463195, 0.00370927, 0.0031137, 0.00261141, 0.00215499, 0.00174491, 0.00138268, 0.00106731, 0.000798828, 0.00057785, 0.00040336, 0.00027161, 0.000176535, 0.00011092, 6.75502e-05, 4.00323e-05, 2.32123e-05, 1.32585e-05, 7.51611e-06, 4.25902e-06, 2.42513e-06, 1.39077e-06, 8.02452e-07, 4.64159e-07, 2.67845e-07, 1.5344e-07, 8.68966e-08, 4.84931e-08, 2.6606e-08, 1.433e-08]
            ),
        ## mix_2015_25ns_HiLum_PoissonOOTPU_cfi https://github.com/cms-sw/cmssw/blob/master/SimGeneral/MixingModule/python/mix_2015_25ns_HiLum_PoissonOOTPU_cfi.py#L24-L77
        "Spring15_25ns": (
            np.linspace(0., 52., 53),
            [4.8551E-07, 1.74806E-06, 3.30868E-06, 1.62972E-05, 4.95667E-05, 0.000606966, 0.003307249, 0.010340741, 0.022852296, 0.041948781, 0.058609363, 0.067475755, 0.072817826, 0.075931405, 0.076782504, 0.076202319, 0.074502547, 0.072355135, 0.069642102, 0.064920999, 0.05725576, 0.047289348, 0.036528446, 0.026376131, 0.017806872, 0.011249422, 0.006643385, 0.003662904, 0.001899681, 0.00095614, 0.00050028, 0.000297353, 0.000208717, 0.000165856, 0.000139974, 0.000120481, 0.000103826, 8.88868E-05, 7.53323E-05, 6.30863E-05, 5.21356E-05, 4.24754E-05, 3.40876E-05, 2.69282E-05, 2.09267E-05, 1.5989E-05, 4.8551E-06, 2.42755E-06, 4.8551E-07, 2.42755E-07, 1.21378E-07, 4.8551E-08]
            )
        }

def getHist(fName, hName="pileup"):
    from bamboo.root import gbl
    tf = gbl.TFile.Open(fName)
    if not tf:
        raise RuntimeError("Could not open file '{0}'".format(fName))
    hist = tf.Get(hName)
    if not hist:
        raise RuntimeError("No histogram with name '{0}' found in file '{1}'".format(hName, fName))
    return tf, hist

def normAndExtract(hist, norm=1.):
    nB = hist.GetNbinsX()
    xAx = hist.GetXaxis()
    if norm:
        hist.Scale(norm/(hist.Integral()*(xAx.GetXmax()-xAx.GetXmin())/nB))
    bEdges = np.array([ xAx.GetBinLowEdge(i) for i in range(1,nB+1) ]+[ xAx.GetBinUpEdge(nB) ])
    contents = np.array([ hist.GetBinContent(i) for i in range(1,nB+1) ])
    return bEdges, contents

def getRatio(numBins, numCont, denBins, denCont):
    ## use numerator for output format
    if not all(db in numBins for db in denBins):
        raise RuntimeError("Numerator (data) needs to have at least the bin edges that are the denominator (MC)")
    ## ratios for the common range
    xMinC, xMaxC = denBins[0], denBins[-1]
    inMn = np.where(numBins == xMinC)[0][0]
    inMx = np.where(numBins == xMaxC)[0][0]
    ratio = np.zeros((inMx-inMn,))
    di = 0
    for ni in range(inMn, inMx):
        if numBins[ni+1] > denBins[di+1]:
            di += 1
        assert ( denBins[di] <= numBins[ni] ) and ( numBins[ni+1] <= denBins[di+1] )
        ratio[ni-inMn] = numCont[ni]/denCont[di]
    bR = np.array(numBins[inMn:inMx+1])
    ## extend range of outside ratio bins until end of numerator ranges
    bR[0] = numBins[0]
    bR[-1] = numBins[-1]
    return bR, ratio

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Produce a BinnedValues-JSON file for pileup reweighting, using data pileup distributions obtained with `pileupCalc.py -i analysis-lumi-json.txt --inputLumiJSON pileup-json.txt --calcMode true --minBiasXsec MBXSECINNB --maxPileupBin NMAX --numPileupBins N outname.root` (see also https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II)")
    parser.add_argument("-o", "--output", default="puweights.json", type=str, help="Output file name")
    parser.add_argument("--mcprofile", default="Moriond17_25ns", help="Pileup profile used to generate the MC sample (use --listmcprofiles to see the list of defined profiles)")
    parser.add_argument("--listmcprofiles", action="store_true", help="list the available MC pileup profiles")
    parser.add_argument("--nominal", type=str, help="File with the data (true) pileup distribution histogram assuming the nominal minimum bias cross-section value")
    parser.add_argument("--up", type=str, help="File with the data (true) pileup distribution histogram assuming the nominal+1sigma minimum bias cross-section value")
    parser.add_argument("--down", type=str, help="File with the data (true) pileup distribution histogram assuming the nominal-1sigma minimum bias cross-section value")
    parser.add_argument("--rebin", type=int, help="Factor to rebin the data histograms by")
    parser.add_argument("--makePlot", action="store_true", help="Make a plot of the PU profiles and weights (requires matplotlib)")
    args = parser.parse_args()
    if args.listmcprofiles:
        print("The known PU profiles are: {0}".format(", ".join(repr(k) for k in mcPUProfiles)))
    else:
        if args.mcprofile not in mcPUProfiles:
            raise ValueError("No MC PU profile with tag '{0}' is known".format(args.mcprofile))
        if not args.nominal:
            raise RuntimeError("No --nominal argument")

        mcPUBins, mcPUVals = mcPUProfiles[args.mcprofile]
        if len(mcPUBins) != len(mcPUVals)+1:
            print(len(mcPUBins), len(mcPUVals))

        fNom, hNom = getHist(args.nominal)
        if args.rebin:
            hNom.Rebin(args.rebin)
        nomBins, nomCont = normAndExtract(hNom)
        ratioBins, nomRatio = getRatio(nomBins, nomCont, mcPUBins, mcPUVals)

        upCont, upRatio, downCont, downRatio = None, None, None, None
        if bool(args.up) != bool(args.down):
            raise ValueError("If either one of --up and --down is specified, both should be")
        if args.up and args.down:
            fUp, hUp = getHist(args.up)
            if args.rebin:
                hUp.Rebin(args.rebin)
            upBins, upCont = normAndExtract(hUp)
            #if not all(ub == nb for ub,nb in zip(upBins, nomBins)):
            #    raise RuntimeError("Up-variation and nominal binning is different: {0} vs {1}".format(upBins, nomBins))
            _, upRatio = getRatio(upBins, upCont, mcPUBins, mcPUVals)
            fDown, hDown = getHist(args.down)
            if args.rebin:
                hDown.Rebin(args.rebin)
            downBins, downCont = normAndExtract(hDown)
            #if not all(db == nb for db,nb in zip(downBins, nomBins)):
            #    raise RuntimeError("Up-variation and nominal binning is different: {0} vs {1}".format(upBins, nomBins))
            _, downRatio = getRatio(downBins, downCont, mcPUBins, mcPUVals)

        out = {
              "dimensin" : 1
            , "variables" : ["NumTrueInteractions"]
            , "binning" : {"x": list(ratioBins)}
            , "error_type" : "absolute"
            , "data" : [
                { "bin" : [ratioBins[i], ratioBins[i+1]]
                , "value" : nomRatio[i]
                , "error_low"  : (nomRatio[i]-downRatio[i] if downRatio is not None else 0.)
                , "error_high" : (upRatio[i]-nomRatio[i] if upRatio is not None else 0.)
                } for i in range(nomRatio.shape[0])
                ]
            }
        with open(args.output, "w") as outF:
            json.dump(out, outF)

        if args.makePlot:
            try:
                from matplotlib import pyplot as plt
            except Exception as ex:
                logger.error("matplotlib could not be imported, so no plot will be produced")
                import sys; sys.exit(0)

            fig,(ax,rax) = plt.subplots(2,1,sharex=True)
            rax.semilogy()
            #rax = ax.twinx()
            dBinCenters = .5*(mcPUBins[:-1]+mcPUBins[1:])
            nBinCenters = .5*(nomBins[:-1]+nomBins[1:])
            rBinCenters = .5*(ratioBins[:-1]+ratioBins[1:])
            ax.hist(dBinCenters, bins=mcPUBins, weights=mcPUVals, histtype="step", label="MC")
            ax.hist(nBinCenters, bins=nomBins, weights=nomCont, histtype="step", label="Nominal", color="k")
            rax.hist(rBinCenters, bins=ratioBins, weights=nomRatio, histtype="step", color="k")
            if upCont is not None:
                ax.hist(nBinCenters, bins=nomBins, weights=upCont, histtype="step", label="Up", color="r")
                ax.hist(nBinCenters, bins=nomBins, weights=downCont, histtype="step", label="Down", color="b")
                rax.hist(rBinCenters, bins=ratioBins, weights=upRatio, histtype="step", color="r")
                rax.hist(rBinCenters, bins=ratioBins, weights=downRatio, histtype="step", color="b")
            rax.axhline(1.)
            ax.legend()
            rax.set_ylim(.1, 2.)
            plt.show()

if __name__ == "__main__":
    main()
