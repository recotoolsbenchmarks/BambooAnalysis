import dataframe.treefunctions as op
from dataframe.treedecorators import decorateTTW
from dataframe.plots import Selection, Plot, EquidistantBinning
from dataframe.dataframebackend import DataframeBackend
import ROOT
import yaml
with open("ttW.yaml") as df:
    ttWdesc = yaml.load(df)
f = ROOT.TFile.Open("ttW_sample.root")
t_ = f.Get("t")
t = decorateTTW(t_, description=ttWdesc)
be, noSel = DataframeBackend.create(t)
