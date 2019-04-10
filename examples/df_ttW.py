import dataframe.treefunctions as op
from dataframe.treedecorators import decorateTTW
from dataframe.plots import Selection, Plot, EquidistantBinning
from dataframe.dataframebackend import DataframeBackend
from cppyy import gbl
import yaml
with open("ttW.yaml") as df:
    ttWdesc = yaml.load(df)
f = gbl.TFile.Open("ttW_sample.root")
t_ = f.Get("t")
t = decorateTTW(t_, description=ttWdesc)
be, noSel = DataframeBackend.create(t)
