# Bamboo: A high-level HEP analysis library for ROOT::RDataFrame

The [`ROOT::RDataFrame`](https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html)
class provides an efficient and flexible way to process per-event information
(stored in a `TTree`) and e.g. aggregate it into histograms.
With the typical pattern of storing object arrays as a structure of arrays
(variable-sized branches with a common prefix in the names and length),
the expressions that are typically needed for a complete analysis quickly become
cumbersome to write, (with indices to match, repeated sub-expressions etc.) -
as an example, imagine the expression needed to calculate the invariant mass of
the two leading muons from a NanoAOD (which stores momenta with pt, eta and phi
branches - one way is to construct LorentzVector objects, sum and evaluate the
invariant mass.
Next imagine doing the same thing with the two highest-PT jets that have a b-tag
and are not within some cone of the two leptons you already selected in another
way (while keeping the code maintainable enough to allow for passing jet momenta
with a systematic variation applied).

This library attempt to solve this problem by automatically constructing
lightweight python wrappers based on the structure of the `TTree`,
which allow to construct such expression with high-level code, similar to the
language that is commonly used to discuss and describe them.
In addition, as much as possible of the technical code needed for processing the
dataset, saving the results etc. is kept in common modules and classes,
such that only high-level declarative-style code remains in the
analysis module, making it easier to read, and compact enough that all relevant
choices and definitions fit in one manageable file.

Bamboo provides tree decorators (branch proxy objects that emulate the
underlying values etc.) and helper methods to build an object representation
of expressions from them, e.g. invariant masses, but also 'the jet with the
highest b-tagger output' or 'the highest-PT opposite-sign same-flavour lepton
pair passing some cuts'. Since such expressions are composable, a lot of ground
can be covered with a few powerful operations.

Building selections, plots etc. with such expressions is analysis-specific, but
the mechanics of loading data samples, processing them locally or on a batch
system (and merging the output of that), combining the outputs for different
samples in an overview etc. is very similar over a broad range of use cases.
Therefore a common implementation of these is provided, which can be used by
extending a base class (to fill histograms and make stacked plots, a class
needs to be written with a method that returns a list of 'plot' objects -
each essentially a combination of an x-axis variable, selection, and weight
to apply to every event - and a configuration file that specifies which datasets
should be processed, and how they should be stacked).

## Installation and setup

Bamboo only depends on python3 (with pip/setuptools to install PyYAML if needed)
and a recent version of ROOT (at least 6.14, for the necessary features of
RDataFrame).
On ingrid and lxplus (or any machine with cvmfs), an easy way to get such
a recent version of ROOT is through a CMSSW release that depends on it (`10_4`
has a recent enough ROOT, but no python3 yet), or from the lcgsoft distribution,
e.g.
```bash
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```
(the second command creates a [virtual environment](https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments)
to install python packages in, after installation it is sufficient to run two
other commands, to pick up the correct base system and then the installed
packages).

Bamboo can be installed in the virtual environment with pip, minimally either
one of
```bash
pip install git+https://cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git
pip install git+ssh://git@cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git
```
but in the current development stage it may be useful to install from
a local clone, such that you can use it to test and propose changes, using
```bash
git clone -o upstream git+ssh://git@cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git /path/to/your/bambooclone
pip install /path/to/your/bambooclone ## e.g. ./bamboo (not bamboo - a package with that name exists)
```
(you will need to upgrade with `pip install --upgrade /path/to/your/bambooclone`
still because installing in editable mode does not work well with including C++
libraries as extensions).

For combining the different histograms in stacks and producing pdf or png files,
which is used in many analyses, the [plotIt](https://github.com/cp3-llbb/plotIt)
tool is used. It can be installed with
```bash
git clone -o upstream https://github.com/cp3-llbb/plotIt.git /path/to/your/plotitclone
cd /path/to/your/plotitclone
cd external
./build-external.sh
cd -
BOOST_ROOT=$CMAKE_PREFIX_PATH make -j4
cp plotIt bamboovenv/bin
```
where the last command copies the `plotIt` executable inside the virtualenv
executable directory such that it is picked up automatically (alternatively,
its path can be passed to `bambooRun` with the `--plotIt` command-line option).

Putting the above commands together, the following should give you a virtual
environment with bamboo, and a clone of bamboo and plotIt in case you need to
modify them (they can be updated with `git pull` and `pip install --upgrade`):
```bash
mkdir bamboodev
cd bamboodev
source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
git clone -o upstream git+ssh://git@cp3-git.irmp.ucl.ac.be/pdavid/bamboo.git
pip install ./bamboo
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
cd plotIt/external
./build-external.sh
cd ..
BOOST_ROOT=$CMAKE_PREFIX_PATH make -j4
cd ..
cp plotIt/plotIt bamboovenv/bin
```

Now you can run a few simple tests on a CMS NanoAOD (on ingrid you could use
`/home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root`)
to see if the installation was successful.
First, we can pretend we are a 'worker' task, which processes trees and outputs
a file with histograms, with a [test module](examples/nanozmumu.py):
```bash
bambooRun -m /path/to/your/clone/examples/nanozmumu.py:NanoZMuMu --distributed=worker /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root -o testh1.root
```
(`--distributed=worker` is needed to interpret the positional arguments as
input file names, in sequential mode (no `--distributed` option) and for
the driver task (`--distributed=driver`) the positional argument is reserved
for a json/yaml file that contains more information, such as input file
locations for several samples, normalisation etc. - there are a few examples).

A more complete example would run from an `analysis.yml` file (copy it to
`bamboo/examples` because [`test_nanozmm1.yml`](examples/test_nanozmm1.yml)
specifies it as a local file with relative path):
```bash
cp /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root bamboo/examples
bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu bamboo/examples/test_nanozmm1.yml --envConfig=bamboo/examples/ingrid.ini -o test_nanozmm1
```
if all went well, you should have a dimuon Z peak plot in
`test_nanozmm1/plots/dimu_M.pdf`. To run on slurm add `--distributed=driver`.

Passing the `--envConfig` option can in practice be avoided by copying the
appropriate file to `~/.config/bamboorc`. It is necessary to pick up the
configuration of the computing environment (files for ingrid and lxplus are
included), e.g. how to access the file storage, which batch submission system
to use (currently slurm and HTCondor are supported). Bamboo tries to find it
first from the `--envConfig` option, then from `bamboo.ini` in the current
directory, then `$XDG_CONFIG_HOME/bamboorc` (which typically resolves to
`~/.config/bamboorc`).

## User guide

This section contains some more information on doing your analysis with bamboo.
It assumes you have successfully installed it following the instructions in the
previous section.

The first thing to make sure is that bamboo can work with your trees.
A first set of decorators for CMS NanoAOD is included in
[`bamboo.treedecorators.decorateNanoAOD`](bamboo/treedecorators.py#L139);
to make stacked histogram plots from them it is sufficient to make your
analysis module inherit from
[`bamboo.analysismodules.NanoAODHistoModule`](bamboo/analysismodules.py#L248)
(which calls this method from its `prepareTrees` method).
Other types of trees can be included in a similar way, but a bit of development
is needed to provided a convenient way to do so (help welcome).

Then, you will need to provide the two really analysis-specific parts: a
configuration file that lists the samples and plot configuration, and a module
specifying what to fill in the histograms. Let's start with the latter.

### Analysis YAML file format

The analysis configuration file should be in the [YAML](https://yaml.org)
format. This was chosen because it can easily be parsed while also being very
readayble (see the [YAML Wikipedia page](https://en.wikipedia.org/wiki/YAML)
for some examples and context) - it essentially becomes a nested dictionary,
which can also contain lists.

Two top-level keys are currently required: `tree` with the name of the `TTree`
inside the file (e.g. `tree: Events` for NanoAOD), `samples`.
For stacked histogram plots, a `plotIt` section should also be specified (the
[`bamboo.analysisutils.runPlotIt`](bamboo/analysismodules.py#L144) method will
insert the `files` and `plots` sections and run plotIt with the resulting
configuration).

Each entry in the `samples` dictionary (the keys are the names of the samples)
is another dictionary. The files to processed can be specified directly as a
list under `files` (with paths relative to the location of the config file,
which is useful for testing), or absolute paths/urls (e.g. xrootd).
If `files` is a string, it is taken as a file with a list of such paths/urls.
For actual analyses, however, samples will be retrieved from a database, e.g.
[DAS](https://cmsweb.cern.ch/das/) or
[SAMADhi](https://cp3.irmp.ucl.ac.be/samadhi/index.php) (support for the latter
still needs to be implemented).
In that case, the database path or query can be specified under `db`, e.g.
`db: das:/SingleMuon/Run2016E-Nano14Dec2018-v1/NANOAOD`.
If both `db` and `files` are specified, and `files` is a string, it is taken as
the path of a cache file with the results from the query: if it does not exist
the query is performed and the result written to the cache file; if it does
exist the list of files is read directly from there. The latter can be
overridden with the `--redodbqueries` option. If in addition the
`--overwritesamplefilelists` option is specified, the results will be saved
(even if the files exist); the cache can also be refreshed by removing the
cache files.

For data, it is usually necessary to specify a json file to filter the good
luminosity blocks (and a run range to consider from it, for efficiency).
If an url is specified for the json file, the file will be downloaded
automatically (and added to the input sandbox for the worker tasks, if needed).

For the formatting of the stack plots, each sample needs to be in a group (e.g.
'data' for data etc.), which will be taken together as one contribution.
An `era` key is also foreseen (to make 2016/2017/2018/combined plots) - but
it is currently ignored.

For the normalization of simulated samples in the stacks, the number of
generated evens and cross-section are also needed. The latter should be
specified as `cross-section` with the sample (in the same units as `luminosity`
in the `configuration` subsection of `plotIt`), the former can be computed from
the input files. For this, the
[`bamboo.analysismodules.HistogramsModule`](bamboo/analysismodules.py#L189)
base class will call the `mergeCounters` method when processing the samples, and
the `readCounters` to read the values from the results file - for NanoAOD
the former merges the `Runs` trees and saves the results, while the latter
performs the sum of the branch with the name specified under `generated-events`.

All together, typical data and MC sample entries would look like
```yaml
  SingleMuon_2016E:
    db: das:/SingleMuon/Run2016E-Nano14Dec2018-v1/NANOAOD
    files: dascache/SingleMuon_2016E.dat
    run_range: [276831, 277420]
    certified_lumi_file: https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt
    era: 2016
    group: data

  DY_high_2017:
    db: das:/DYJetsToLL_M-50_TuneCP5_13TeV-amcatnloFXFX-pythia8/RunIIFall17NanoAODv4-PU2017_12Apr2018_Nano14Dec2018_new_pmx_102X_mc2017_realistic_v6_ext1-v1/NANOAODSIM
    files: dascache/DY_high_2017.dat
    era: 2017
    group: DY
    cross-section: 5765.4
    generated-events: genEventSumw
```

### Analysis module

For an analysis module to be run with `bambooRun`, it in principle only needs
a constructor that takes an argument with command-line arguments, and a `run`
method. [`bamboo.analysismodules`](bamboo/analysismodules.py) provides a more
interesting base class [`bamboo.AnalysisModule`](bamboo/analysismodules.py#L39)
that provides a lot of common functionality (most notably: parsing the analysis
configuration, running sequentially or distributed (and also as worker task in
the latter case), and provides `addArgs(parser)`, `initialize()`,
`processTrees(inputFiles, outputFile, tree=None, certifiedLumiFile=None, runRange=None)`,
`postProcess(taskList, config=None, workdir=None, resultsdir=None)` and
`interact()` interface member methods that should be further specified by
subclasses.

[`bamboo.analsysimodules.HistogramsModule`](bamboo/analysismodules.py#L189)
does this for the stacked histogram plots, composing `processTrees` from
`prepareTree(tree)` and `definePlots(decoratedTree, systVar="nominal")`,
while taking the JSON lumi block mask and counter merging into account.
It also calls the `plotIt` calls from `postProcess` (with the plots list,
has all required information).
[`bamboo.analysismodules.NanoAODHistoModule`](bamboo/analysismodules.py#L248)
supplements this with the decorations and counter merging and reading for
NanoAOD, such that all the final module needs to do is defining plots and
selections, as in [this example](example/nanozmumu.py).
This layered structure is used such that code can be maximally reused for other
types of trees.

For the code inside the module, the example is also very instructive:
```python
    def definePlots(self, t, noSel, systVar="nominal"):
        from bamboo.plots import Plot, EquidistantBinning
        import bamboo.treefunctions as op

        plots = []

        twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(t.Muon) > 1 ])
        plots.append(Plot.make1D("dimu_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSel,
                EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

        return plots
```
The key classes are defined in [`bamboo.plots`](bamboo/plots.py):
[`Plot`](bamboo/plots.py#L55) and [`Selection`](bamboo/plots.py#L115).
The latter represents a consistent set of selection requirements (cuts) and
weight factors (e.g. to apply corrections). Selections are defined by refining
a "root selection" with additional cuts and weights, and each should have a
unique name (an exception is raised at construction otherwise).
The root selection allows to do some customisation upfront, e.g. the applying
the JSON luminosity block mask for data. A plot object refers to a selection,
and specifies which variable(s) to plot, with which binning(s), labels, options
etc. (the `plotOpts` dictionary is copied directly into the plot section of the
plotIt configuration file).

## Specifying cuts, weight, and variables: expressions

The first argument to the `definePlots` method is the "decorated" tree - a
proxy object from which expressions can be derived. Sticking with the NanoAOD
example, `t.Muon` is another proxy object for the muon collection (similarly
for the other objects), `t.Muon[0]` retrieves the leading-PT muon proxy, and
`t.Muon[0].p4` its momentum fourvector.
The proxies are designed to behave as much as possible as the value types they
correspond to (you can get an item from a list, an attribute from an object,
you can also work with numerical values, e.g.
`t.Muon[0].p4.Px()+t.Muon[1].p4.Px()`) but for some more complex operations,
specific functions are needed. These are as much as possible defined in the
[`bamboo.treefunctions`](python/treefunctions.py) module, so it is recommended
to
```python
import bamboo.treefunctions as op
```
(or any other name, as you prefer) in your analysis module.

Ideally, the decorated tree and the `bamboo.treefunctions` are all you ever need
to import and know about the decorations. Therefore the best way to proceed now
is get a decorated tree inside an IPython shell and play around. For 
[`bamboo.analysismodules.HistogramsModule`](bamboo/analysismodules.py#L189) this
can always be done by passing the `--interactive` flag, with either one of
(depending on if you copied the NanoAOD test file above)
```bash
bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive --distributed=worker /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root
bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive bamboo/examples/test_nanozmm.yml [ --envConfig=bamboo/examples/ingrid.ini ] -o int1
```
The decorated tree is in the `tree` variable (the original `TChain` is in `tup`)
and the `bamboo.treefunctions` module is there as `op` (the `c_...` methods
construct a constant, whereas the `rng_...` methods work on a collection and
return a single value (note: all these methods should be documented), whereas
the `op.select` method returns a reduced collection (internally, only a list of
indices to the passing objects is created, and the result is a proxy that uses
this list). Some of the `rng_...` methods are extremely powerful, e.g.
`rng_find` and `rng_max_element_by`.
The proxy classes are generated on-the-fly with all branches as attributes, so
tab-completion can be used to have a look at what's there:
```python
In [1]: tree.<TAB>
  tree.CaloMET                           tree.SoftActivityJetHT10
  tree.Electron                          tree.SoftActivityJetHT2
  tree.FatJet                            tree.SoftActivityJetHT5
  tree.Flag                              tree.SoftActivityJetNjets10
  tree.HLT                               tree.SoftActivityJetNjets2
  tree.HLTriggerFinalPath                tree.SoftActivityJetNjets5
  tree.HLTriggerFirstPath                tree.SubJet
  tree.Jet                               tree.Tau
  tree.L1Reco_step                       tree.TkMET
  tree.MET                               tree.TrigObj
  tree.Muon                              tree.deps
  tree.OtherPV                           tree.event
  tree.PV                                tree.fixedGridRhoFastjetAll
  tree.Photon                            tree.fixedGridRhoFastjetCentralCalo
  tree.PuppiMET                          tree.fixedGridRhoFastjetCentralNeutral
  tree.RawMET                            tree.luminosityBlock
  tree.SV                                tree.op
  tree.SoftActivityJet                   tree.run
  tree.SoftActivityJetHT                                                        

In [1]: anElectron = tree.Electron[0]

In [2]: anElectron.<TAB>
   anElectron.charge                   anElectron.eInvMinusPInv            anElectron.mvaSpring16HZZ_WPL
   anElectron.cleanmask                anElectron.energyErr                anElectron.mvaTTH
   anElectron.convVeto                 anElectron.eta                      anElectron.op
   anElectron.cutBased                 anElectron.hoe                      anElectron.p4
   anElectron.cutBased_HEEP            anElectron.ip3d                     anElectron.pdgId
   anElectron.cutBased_HLTPreSel       anElectron.isPFcand                 anElectron.pfRelIso03_all
   anElectron.deltaEtaSC               anElectron.jet                      anElectron.pfRelIso03_chg
   anElectron.dr03EcalRecHitSumEt      anElectron.lostHits                 anElectron.phi
   anElectron.dr03HcalDepth1TowerSumEt anElectron.mass                     anElectron.photon
   anElectron.dr03TkSumPt              anElectron.miniPFRelIso_all         anElectron.pt
   anElectron.dxy                      anElectron.miniPFRelIso_chg         anElectron.r9
   anElectron.dxyErr                   anElectron.mvaSpring16GP            anElectron.sieie
   anElectron.dz                       anElectron.mvaSpring16GP_WP80       anElectron.sip3d
   anElectron.dzErr                    anElectron.mvaSpring16GP_WP90       anElectron.tightCharge
   anElectron.eCorr                    anElectron.mvaSpring16HZZ           anElectron.vidNestedWPBitmap
```
