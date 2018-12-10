# Bamboo: A high-level HEP analysis library for ROOT::RDataFrame

The [`ROOT::RDataFrame`](https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html)
class provides an efficient and flexible way to process per-event information
(stored in a `TTree`) and e.g. aggregate it into histograms.
With the typical pattern of storing object arrays as a structure of arrays
(variable-sized branches with a common prefix in the names), the expressions
that are typically needed for a complete analysis quickly become cumbersome
to write, (with indices to match, repeated sub-expressions etc.).

This library attempt to alleviate this problem by automatically constructing
lightweight python wrappers based on the structure of the `TTree`,
such that selections, plotted variables etc. can be written as
high-level expressions that are much closer to the language commonly used to
discuss and describe them.
Making the analysis code more compact allows to keep all relevant choices and
definitions close together, which helps in keeping the overview.

As an example, consider a CMS NanoAOD tree, where particle momenta are stored
split into `X_pt[]`, `X_eta[]`, `X_phi[]` and `X_mass[]` branches.
Calculating the invariant mass of two particles would entail something like
`(LV(X_pt[i1], X_eta[i1], X_phi[i1], X_mass[i1]) +
LV(X_pt[i2], X_eta[i2], X_phi[i2], X_mass[i2])).M()` without this library,
but could be simplified to `(X[i1].p4+X[i2].p4).M()` with it, by defining a
trivial `p4` property to construct this on the fly.

More complex expressions can also be constructed by building expressions using
the columns of the `TTree`: the proxy instances return an object
representation of the expression instead of the actual result of a calculation,
but otherwise mimic the number, array or object they correspond to.

A particularly powerful class of operations are those that work on columns that
contain arrays, to extract for instance the sum/average, the first/best item
according to some criterion.

``` python
import bamboo.treefunctions as op
## t is the decorated NanoAOD
mostCentralMuon = op.rng_min_element_by(t.Muon, lambda mu : op.abs(mu.eta))
## the mostCentralMuon proxy object can be used to draw e.g. its pt
mostCentralMuon.p4.Pt()
```

The `plots` module also defines `Plot` and `Selection` classes, with the
latter representing a consistent set of selection requirements (cuts) and
weight factors (e.g. to apply corrections), and the former such a `Selection`,
one or more variables to use for the x-axis (y-axis, ...), the corresponding
binnings, and other information to make up the graph (axis titles, line/area
styles etc.). These are intended to be the basic building blocks of analysis
workflows.

## Getting started

Bamboo needs only python3 and a recent version of ROOT (at least 6.14,
for the necessary features of RDataFrame). On ingrid and lxplus (or any machine
with cvmfs), an easy way to get such a recent version of ROOT is through
a CMSSW release that depends on it (`10_4` or later), or from the lcgsoft
distribution, e.g.
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```

The package can be installed with pip, minimally either one of
```bash
% pip install git+https://cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git
% pip install git+ssh://git@cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git
```
but in the current development stage it may be useful to install from
a local clone, such that you can use it to test and propose changes, using
```bash
% git clone -o upstream git+ssh://git@cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git /path/to/your/clone
% pip install /path/to/your/clone
```
(you will need to upgrade with `pip install --upgrade /path/to/your/clone` still
because installing in editable mode does not work well with extensions etc.).

On ingrid-ui1, you can run the following test to check the installation:
```bash
% bambooRun -m /home/ucl/cp3/pdavid/bambootest/nanozmumu.py:NanoZMuMu --worker /home/ucl/cp3/pdavid/bambootest/NanoAOD_SingleMu_test.root
```