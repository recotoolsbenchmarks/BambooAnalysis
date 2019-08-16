User guide
==========

This section contains some more information on doing your analysis with bamboo.
It assumes you have successfully installed it following the instructions in the
:doc:`previous section<install>`.

The first thing to make sure is that bamboo can work with your trees.
With CMS NanoAOD, many analysis use the same (or a very similar) tree format,
which is why a set of decorators for is included in
:py:func:`bamboo.treedecorators.decorateNanoAOD`; to make stacked histogram
plots from them it is sufficient to make your analysis module inherit from
:py:class:`bamboo.analysismodules.NanoAODHistoModule`
(which calls this method from its
:py:meth:`~bamboo.analysismodules.NanoAODHistoModule.prepareTrees` method).
Other types of trees can be included in a similar way, but a bit of development
is needed to provided a more convenient way to do so (help welcome).

Running bambooRun
-----------------

The ``bambooRun`` executable script can be used to run over some samples and
derive other samples or histograms from them. It needs at least two arguments: a
:ref:`python module<uganalysismodule>`, which will tell it what to do for each event,
and a :ref:`configuration file<uganalysisyaml>` with a list of samples to process,
plot settings etc.

Typically, ``bambooRun`` would be invoked with

.. code-block:: sh

   bambooRun -m module-specification config-file-path -o my-output

where ``module-specification`` is of the format ``modulename:classname``, and
``modulename`` can be either a file path like ``somedir/mymodule.py`` or an
importable module name  like ``myanalysispackage.mymodule``.
This will construct an instance of the specified module, passing it any
command-line arguments that are not used directly by ``bambooRun``, and run it.

The default base module (:py:class:`bamboo.analysismodules.AnalysisModule`, see
:ref:`below<uganalysismodule>`) provides a number of convenient command-line
options (and individual modules can add more by implementing the
:py:meth:`~bamboo.analysismodules.AnalysisModule.addArgs` method).

* the ``-h`` (``--help``) switch prints the complete list of supported
  options and arguments, including those defined by the module if used with
  ``bambooRun -h -m mymodule``
* the ``-o`` (``--output``) option can be used to specify the base directory for
  output files
* the ``-v`` (``--verbose``) switch will produce more output messages, and also
  print the full C++ code definitions that are passed to RDataFrame_ (which is
  very useful for debugging)
* the ``-i`` (``--interactive``) switch will only load one file and launch an
  IPython terminal, where you can have a look at its structure and test
  expressions
* the ``--maxFiles`` option can be used to specify a maximum number of files
  to process for each sample, e.g. ``--maxFiles=1`` to check that the module
  runs correctly in all cases before submitting to a batch system

The usual mode of operation is to parse the analysis configuration file,
execute some code for every entry in each of the samples, and then perform some
actions on the aggregated results (e.g. draw histograms).
Since the second step is by far the most time-consuming, but can be performed
indepently for different samples (and even entries), it is modeled as a list of
tasks (which may be run in parallel), after which a postprocessing step takes
the results ad combines them (this can also be run separately, using
the results of previously run tasks, assuming these did not change).

More concretely, for histogram stack plots the tasks produce histograms while
the postprocessing step runs plotIt_, so with the ``--onlypost`` option the
normalization, colors, labels etc. can be changed without reprocessing the
samples.
Passing the ``--distributed=driver`` option will submit the independent tasks to
a batch scheduler (currently HTCondor and Slurm are supported) instead of
running them sequentially, wait for the results to be ready, and combine them
(the worker tasks will run the same module, but with ``--distributed=worker``
and the actual input and results file names as input and output arguments).
By default one batch job is submitted for each input sample, unless there is
a ``split`` entry different from one for the sample, see
:ref:`below<uganalysisyaml>` for the precise meaning.

.. _ugenvconfig:

Computing environment configuration file
''''''''''''''''''''''''''''''''''''''''

For some features such as automatically converting logical filenames from DAS
to physical filenames at your local T2 storage (or falling back to xrootd),
submitting to a batch cluster etc., some information about the computing
resources and environment is needed.
In order to avoid proliferating the command-line interface of ``bambooRun``,
these pieces of information are bundled in a file that can be passed in one go
through the ``--envConfig`` option.
If not specified, Bamboo_ will try to read ``bamboo.ini`` in the current
directory, and then ``$XDG_CONFIG_HOME/bamboorc`` (which typically resolves to
``~/.config/bamboorc``).
Since these settings are not expected to change often or much, it is advised to
copy the closest example (e.g. ``examples/ingrid.ini`` or
``examples/lxplus.ini``) to ``~/.config/bamboorc`` and edit if necessary.


.. _uganalysisyaml:

Analysis YAML file format
-------------------------

The analysis configuration file should be in the YAML_ format.  This was chosen
because it can easily be parsed while also being very readable (see the
`YAML Wikipedia page`_ for some examples and context) - it essentially becomes
a nested dictionary, which can also contain lists.

Two top-level keys are currently required: ``tree`` with the name of the TTree_
inside the file (e.g. ``tree: Events`` for NanoAOD), ``samples``.
For stacked histogram plots, a ``plotIt`` section should also be specified (the
:py:func:`bamboo.analysisutils.runPlotIt` method will insert the ``files`` and
``plots`` sections and run plotIt_ with the resulting configuration).

Each entry in the ``samples`` dictionary (the keys are the names of the samples)
is another dictionary. The files to processed can be specified directly as a
list under ``files`` (with paths relative to the location of the config file,
which is useful for testing), or absolute paths/urls (e.g. xrootd).
If ``files`` is a string, it is taken as a file with a list of such paths/urls.
For actual analyses, however, samples will be retrieved from a database, e.g.
DAS_ or SAMADhi_ (support for the latter still needs to be implemented).
In that case, the database path or query can be specified under ``db``, e.g.
``db: das:/SingleMuon/Run2016E-Nano14Dec2018-v1/NANOAOD``.
If both ``db`` and ``files`` are specified, and ``files`` is a string, it is
taken as the path of a cache file with the results from the query: if it does
not exist the query is performed and the result written to the cache file; if it
does exist the list of files is read directly from there. The latter can be
overridden with the ``--redodbqueries`` option. If in addition the
``--overwritesamplefilelists`` option is specified, the results will be saved
(even if the files exist); the cache can also be refreshed by removing the
cache files.

.. tip::

   Samples in DAS and SAMADhi rarely change, and reading a local file is almost
   always faster than doing queries (and does not require a grid proxy etc.),
   so especially when using many samples from these databases it is recommended
   to cache the file lists resulting from these results, by specifying a path
   under ``files`` as well as a database path under ``db``
   (see below for an example).

For data, it is usually necessary to specify a json file to filter the good
luminosity blocks (and a run range to consider from it, for efficiency).
If an url is specified for the json file, the file will be downloaded
automatically (and added to the input sandbox for the worker tasks, if needed).

For the formatting of the stack plots, each sample needs to be in a group (e.g.
'data' for data etc.), which will be taken together as one contribution.
An ``era`` key is also foreseen (to make 2016/2017/2018/combined plots) - but
it is currently ignored.

For the normalization of simulated samples in the stacks, the number of
generated evens and cross-section are also needed. The latter should be
specified as ``cross-section`` with the sample (in the same units as
``luminosity`` in the ``configuration`` subsection of ``plotIt``), the former
can be computed from the input files. For this, the
:py:class:`bamboo.analysismodules.HistogramsModule` base class will call the
``mergeCounters`` method when processing the samples, and the ``readCounters``
method to read the values from the results file - for NanoAOD the former merges
the `Runs` trees and saves the results, while the latter performs the sum of
the branch with the name specified under ``generated-events``.

For large samples, a ``split`` property can be specified, such that the input
files are spread out over different batch jobs.
A positive number is taken as the number of jobs to divide the inputs over,
while a negative number gives the number of files per job (e.g. ``split: 3``
will create three jobs to process the sample, while ``split: -3`` will result
in jobs that process three files each).

All together, typical data and MC sample entries would look like

.. code-block:: yaml

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
       split: 3


.. tip::
    
    It is possible to insert the content of a configuration file into another, e.g. to separate or reuse the plot- and samples-related setings: simply use the syntax ``!include file.yml`` in the exact place where you would like to insert the content of ``file.yml``.


.. _uganalysismodule:

Analysis module
---------------

For an analysis module to be run with ``bambooRun``, it in principle only needs
a constructor that takes an argument with command-line arguments, and a ``run``
method. :py:mod:`bamboo.analysismodules` provides a more interesting base class
:py:class:`~bamboo.analysismodules.AnalysisModule` that provides a lot of common
functionality (most notably: parsing the analysis configuration, running
sequentially or distributed (and also as worker task in the latter case), and
provides :py:meth:`~bamboo.analysismodules.AnalysisModule.addArgs`,
:py:meth:`~bamboo.analysismodules.AnalysisModule.initialize`, 
:py:meth:`~bamboo.analysismodules.AnalysisModule.processTrees`,
:py:meth:`~bamboo.analysismodules.AnalysisModule.postProcess`, and
:py:meth:`~bamboo.analysismodules.AnalysisModule.interact`, interface member
methods that should be further specified by subclasses (see the
:doc:`reference documentation<apiref>` for more details).

:py:class:`~bamboo.analysismodules.HistogramsModule` does this for the
stacked histogram plots, composing
:py:meth:`~bamboo.analysismodules.HistogramsModule.processTrees` from
:py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree` and
:py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots`, while taking
the JSON lumi block mask and counter merging into account.
It also calls the `plotIt` executable from 
:py:meth:`~bamboo.analysismodules.HistogramsModule.postProcess` (with the plots
list and analysis configuration file, it has all required information for that).
:py:class:`~bamboo.analysismodules.NanoAODHistoModule` supplements this with
the decorations and counter merging and reading for NanoAOD,
such that all the final module needs to do is defining plots and selections,
as in the example :py:mod:`examples.nanozmumu`.
This layered structure is used such that code can be maximally reused for other
types of trees.

For the code inside the module, the example is also very instructive:

.. code-block:: python

       def definePlots(self, t, noSel, era=None, sample=None):
           from bamboo.plots import Plot, EquidistantBinning
           from bamboo import treefunctions as op

           plots = []

           twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(t.Muon) > 1 ])
           plots.append(Plot.make1D("dimu_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSel,
                   EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

           return plots

The key classes are defined in :py:mod:`bamboo.plots`:
:py:class:`~bamboo.plots.Plot` and :py:class:`~bamboo.plots.Selection`
(see the :doc:`reference documentation<apiref>` for details).
The latter represents a consistent set of selection requirements (cuts) and
weight factors (e.g. to apply corrections). Selections are defined by refining
a "root selection" with additional cuts and weights, and each should have a
unique name (an exception is raised at construction otherwise).
The root selection allows to do some customisation upfront, e.g. the applying
the JSON luminosity block mask for data.
A plot object refers to a selection, and specifies which variable(s) to plot,
with which binning(s), labels, options etc. (the ``plotOpts`` dictionary is
copied directly into the plot section of the plotIt configuration file).

Histograms corresponding to systematic variations (of scalefactors, collections
etc. |---| see below) are by default generated automatically alongside the
nominal one.
This can however easily be disabled at the level of a
:py:class:`~bamboo.plots.Selection` (and, consequently, all
:py:class:`~bamboo.plots.Selection` instances deriving from it, and all
:py:class:`~bamboo.plots.Plot` instances using it) or a single plot, by passing
``autoSyst=False`` to the :py:func:`~bamboo.plots.Selection.refine` or
:py:func:`~bamboo.plots.Plot.make1D` (or related) method, respectively,
when constructing them; so setting ``noSel.autoSyst = False`` right after
retrieving the decorated tree and root selection would turn disable all
automatic systematic variations.

.. _ugexpressions:

Specifying cuts, weight, and variables: expressions
---------------------------------------------------

The first argument to the
:py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots`
method is the "decorated" tree |---| a proxy object from which expressions
can be derived. Sticking with the NanoAOD example, ``t.Muon`` is another
proxy object for the muon collection (similarly for the other objects),
``t.Muon[0]`` retrieves the leading-pt muon proxy, and ``t.Muon[0].p4``
its momentum fourvector.
The proxies are designed to behave as much as possible as the value types they
correspond to (you can get an item from a list, an attribute from an object,
you can also work with numerical values, e.g.
``t.Muon[0].p4.Px()+t.Muon[1].p4.Px()``) but for some more complex operations,
specific functions are needed. These are as much as possible defined in the
:py:mod:`bamboo.treefunctions` module, see :doc:`Building expressions<treefunctions>`
for an overview of all the available methods.

Ideally, the decorated tree and the :py:mod:`bamboo.treefunctions` module
are all you ever need to import and know about the decorations.
Therefore the best way to proceed now is get a decorated tree
inside an IPython shell and play around.
For :py:mod:`bamboo.analysismodules.HistogramsModule` this can always be done
by passing the ``--interactive`` flag, with either one of
(depending on if you copied the NanoAOD test file above)

.. code-block:: sh

   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive --distributed=worker bamboo/tests/data/DY_M50_2016.root
   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive bamboo/examples/test_nanozmm.yml [ --envConfig=bamboo/examples/ingrid.ini ] -o int1

The decorated tree is in the ``tree`` variable (the original ``TChain`` is in
``tup``) and the :py:mod:`bamboo.treefunctions` module is there as `op`
(the ``c_...`` methods construct a constant, whereas the ``rng_...`` methods
work on a collection and return a single value,
whereas the :py:func:`~bamboo.treefunctions.select` method returns
a reduced collection (internally, only a list of indices to the passing objects
is created, and the result is a proxy that uses this list).
Some of the ``rng_...`` methods are extremely powerful, e.g.
:py:func:`~bamboo.treefunctions.rng_find` and 
:py:func:`~bamboo.treefunctions.rng_max_element_by`.

The proxy classes are generated on the fly with all branches as attributes, so
tab-completion can be used to have a look at what's there:

.. code-block:: python

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

For NanoAOD the content of the branches is documented in the various branches of
`this directory <https://cms-nanoaod-integration.web.cern.ch/integration/>`_,
e.g. NanoAODv4 for
`2016 MC <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc94X2016_doc.html>`_,
`2017 MC <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc94Xv2_doc.html>`_,
`2018 MC <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html>`_, and for
`2016 data <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data94X2016_doc.html>`_,
`2017 data <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data94Xv2_doc.html>`_, and
`2018 data <https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/data101X_doc.html>`_.


.. _ugcutordering:

Ordering selections and plots efficiently
'''''''''''''''''''''''''''''''''''''''''

Internally, Bamboo uses the RDataFrame_ class to process the input samples and
produce histograms or skimmed trees |---| in fact no python code is run while
looping over the events: Bamboo builds up a computation graph when
:py:class:`~bamboo.plots.Selection` and :py:class:`~bamboo.plots.Plot`
objects are defined by the analysis module's
:py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots` method,
RDataFrame_ compiles the expressions for the cuts and variables, and the input
files and events are only looped over once, when the histograms are retrieved
and stored.

In practice, however, there are not only ``Filter``
(:py:class:`~bamboo.plots.Selection`) and ``Fill``
(:py:class:`~bamboo.plots.Plot`) nodes in the computation graph, but also
``Define`` nodes that calculate a quantity based on other columns and make
the result available for downstream nodes to use directly.
This is computationally more efficient if the calculation is complex enough.
Bamboo tries to make a good guess at which (sub-)expressions are worth
pre-calculating (typically operations that require looping over a collection),
but the order in which plots and selections are defined may still help to avoid
inserting the same operation twice in the computation graph.

The main feature to be aware of is that RDataFrame_ makes a node in the
computation graph for every ``Define`` operation, and the defined column can
only be used from nodes downstream of that.
Logically, however, all defined columns needed for plots or sub-selections of
one selection will need to be evaluated for all events passing this selection,
and the most efficient is to do this only once, so ideally all definitions
should be inserted right after the ``Filter`` operation of the selection, and
before any of the ``Fill`` and subsequent ``Filter`` nodes.
This is a bit of a simplification: it is possible to imagine cases where it can
be better to have a column only defined for the sub-nodes that actually use it,
but then it is hard to know in all possible cases where exactly to insert the
definitions, so the current implementation opts for a simple and predictable
solution: on-demand definitions of subexpressions are done when
:py:class:`~bamboo.plots.Plot` and :py:class:`~bamboo.plots.Selection` objects
are constructed, and they update the computation graph node that other nodes
that derive from the same selection will be based on.
A direct consequence of this is that it is usually more efficient to first
define plots for a stage of the selection, and only then define refined
selections based on it |---| otherwise the subselection will be based on the
node without the columns defined for the plots and, in the common case where
the same plots are made at different stages of the selection, recreate nodes
with the same definitions in its branch of the graph.
As an illustration, the pseudocode equivalent of these two cases is

.. code-block:: python

   ## define first subselection then plots
   ## some_calculation(other_columns) is done twice
   if selectionA:
       if selectionB:
          myColumn1 = some_calculation(other_columns)
          myPlot1B = makePlot(myColumn1)
       myColumn2 = some_calculation(other_columns)
       myPlot1A = makePlot(myColumn2)

   ## define first plots then subselection
   ## some_calculation(other_columns) is only done once
   if selectionA:
       myColumn1 = some_calculation(other_columns)
       myPlot1A = makePlot(myColumn1)
       if selectionB:
          myPlot1B = makePlot(myColumn1)

This is why it is advisable to reserve the
:py:meth:`~bamboo.analysismodules.HistogramsModule.definePlots` method of the
analysis module for defining event and object container selections, and define
helper methods that declare the plots for a given selection |---| with a
parameter that is inserted in the plot name to make sure they are unique, if
used to define the same plots for different selection stages, e.g.

.. code-block:: python

   def makeDileptonPlots(self, sel, leptons, uname):
       from bamboo.plots import Plot, EquidistantBinning
       from bamboo import treefunctions as op
       plots = [
            Plot.make1D("{0}_llM".format(uname),
               op.invariant_mass(leptons[0].p4, leptons[1].p4), sel,
               EquidistantBinning(100, 20., 120.),
               title="Dilepton invariant mass",
               plotopts={"show-overflow":False}
               )
       ]
       return plots

   def definePlots(self, t, noSel, era=None, sample=None):
       from bamboo import treefunctions as op

       plots = []

       muons = op.select(t.Muon, lambda mu : op.AND(mu.p4.Pt() > 20., op.abs(mu.p4.Eta() < 2.4)))

       twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(muons) > 1 ])

       plots += self.makeDileptonPlots(twoMuSel, muons, "DiMu")

       jets = op.select(t.Jet, lambda j : j.p4.Pt() > 30.)

       twoMuTwoJetSel = twoMuSel.refine("twoMuonsTwoJets", cut=[ op.rng_len(jets) > 1 ])

       plots += self.makeDileptonPlots(twoMuTwoJetSel, muons, "DiMu2j")

       return plots

Finally, there are some cases where the safest is to force the inclusion of a
calculation at a certain stage, for instance when performing expensive function
calls, since the default strategy is not to precalculate these because there are
many more function calls that are not costly.
A prime example of this is the calculation of modified jet collections with e.g.
an alternative JEC aplied, which is done in a separate C++ module (see below),
and is probably the slowest operation in most analysis tasks.
The definition can be added explicitly under a selection by calling the
:py:meth:`bamboo.analysisutils.forceDefine` method, e.g. with
``forceDefine(t._Jet.calcProd, mySelection)``.


.. _ugrecipes:

Recipes for common tasks
------------------------

.. _ugrecipescalefactors:

Using scalefactors
''''''''''''''''''

Scalefactors |---| CMS jargon for efficiency corrections for MC, typically
binned in lepton or jet kinematic variables |---| can be generalized to
functions that take some properties of a physics object and return a single
floating-point number.
The :py:mod:`bamboo.scalefactors` module provides support for constructing
such callable objects from the JSON format used in the `CP3-llbb framework`_,
see some examples
`here <https://github.com/cp3-llbb/Framework/tree/CMSSW_8_0_6p/data/ScaleFactors>`_
(these JSON files are produced from the txt or ROOT files provided by the POGs
using simple python
`scripts <https://github.com/cp3-llbb/Framework/tree/CMSSW_8_0_6p/scripts>`_).
Like their inputs, the JSON files contain the nominal scale factor as well as
its up and down systematic variations, so the
:py:class:`~bamboo.scalefactors.ScaleFactor` behaves as a callable that takes
a physics object and an optional `variation` keyword argument (technically,
it wraps a C++ object that gets the correct value from the JSON file).

The :py:meth:`~bamboo.scalefactors.get_scalefactor` method constructs such
objects from a nested dictionary such as the one in :py:mod:`bamboo.llbbSF`:
the first key is a tag (as an example: "electron_2015_76", for electrons in
2015 data, analysed with a ``CMSSW_7_6_X`` release) and the second key is an
identifier of the selection they correspond to (e.g. ``id_loose``).
The value inside this dictionary can be either a single path to a JSON file,
or a list of ``(period, path)`` pairs, in case scalefactors for different
running periods need to be combined (the ``periods`` keyword argument to
:py:meth:`~bamboo.scalefactors.get_scalefactor` can be used to select only
a certain set of these periods).
The combination is done by either weighting or randomly sampling from the
different periods, according to the fraction of the integrated luminosity in
each (by passing ``combine="weight"`` or ``combine="sample"``, respectively).
Jet flavour tagging and dilepton (e.g. trigger) scalefactors are can also be
specified by passing tuples of the light, c-jet and b-jet scalefactor paths,
and tuples of first-if-leading, first-if-subleading, second-if-leading,
and second-if-subleading (to be reviewed for NanoAOD) scalefactor paths,
respectively, instead of a single path.

Histogram variations representing the shape systematic uncertainty due to an
uncertainty on the scalefactor values can be automatically produced by passing
a name to the ``systName`` keyword argument of the
:py:meth:`~bamboo.scalefactors.get_scalefactor` method.

As an example, some basic lepton ID and jet tagging scalefactors could be
included in an analysis on NanoAOD by defining

.. code-block:: python

 import bamboo.scalefactors
 from itertools import chain
 import os.path

 # scalefactor JSON files are in ScaleFactors/<era>/ alongside the module
 def localize_myanalysis(aPath, era="2016legacy"):
     return os.path.join(os.path.dirname(os.path.abspath(__file__)), "ScaleFactors", era, aPath)

 # nested dictionary with path names of scalefactor JSON files
 # { tag : { selection : absole-json-path } }
 myScalefactors = {
     "electron_2016_94" : dict((k,localize_myanalysis(v)) for k, v in dict(
         ("id_{wp}".format(wp=wp.lower()), ("Electron_EGamma_SF2D_{wp}.json".format(wp=wp)))
         for wp in ("Loose", "Medium", "Tight") ).items()),
     "btag_2016_94" : dict((k, (tuple(localize_myanalysis(fv) for fv in v))) for k,v in dict(
         ( "{algo}_{wp}".format(algo=algo, wp=wp),
           tuple("BTagging_{wp}_{flav}_{calib}_{algo}.json".format(wp=wp, flav=flav, calib=calib, algo=algo)
               for (flav, calib) in (("lightjets", "incl"), ("cjets", "comb"), ("bjets","comb")))
         ) for wp in ("loose", "medium", "tight") for algo in ("DeepCSV", "DeepJet") ).items())
     }

 # fill in some defaults: myScalefactors and bamboo.scalefactors.binningVariables_nano
 def get_scalefactor(objType, key, periods=None, combine=None, additionalVariables=None, systName=None):
     return bamboo.scalefactors.get_scalefactor(objType, key, periods=periods, combine=combine,
         additionalVariables=(additionalVariables if additionalVariables else dict()),
         sfLib=myScalefactors, paramDefs=bamboo.scalefactors.binningVariables_nano, systName=systName)

and adding the weights to the appropriate :py:class:`~bamboo.plots.Selection`
instances with

.. code-block:: python

 electrons = op.select(t.Electron, lambda ele : op.AND(ele.cutBased >= 2, ele.p4.Pt() > 20., op.abs(ele.p4.Eta()) < 2.5))
 elLooseIDSF = get_scalefactor("lepton", ("electron_2016_94", "id_loose"), systName="elID")
 hasTwoEl = noSel.refine("hasTwoEl", cut=[ op.rng_len(electrons) > 1 ],
               weight=[ elLooseIDSF(electrons[0]), elLooseIDSF(electrons[1]) ])

 jets = op.select(t.Jet, lambda j : j.p4.Pt() > 30.)
 bJets = op.select(jets, lambda j : j.btagDeepFlavB > 0.2217) ## DeepFlavour loose b-tag working point
 deepFlavB_discriVar = { "BTagDiscri": lambda j : j.btagDeepFlavB }
 deepBLooseSF = get_scalefactor("jet", ("btag_2016_94", "DeepJet_loose"), additionalVariables=deepFlavB_discriVar, systName="bTag")
 hasTwoElTwoB = hasTwoEl.refine("hasTwoElTwoB", cut=[ op.rng_len(bJets) > 1 ],
                  weight=[ deepBLooseSF(bJets[0]), deepBLooseSF(bJets[1]) ])

.. _ugrecipepureweighting:

Pileup reweighting
''''''''''''''''''

Pileup reweighting to make the pileup distribution in simulation match the one
in data is very similar to applying a scalefactor, except that the efficiency
correction is for the whole event or per-object |---| so the same code can be
used.
The ``makePUReWeightJSON.py`` script included in bamboo can be used to make
a JSON file with weights out of a data pileup profile obtainedby running
``pileupcalc.py``
(inside CMSSW, see the `pileupcalc documentation`_ for details), e.g. with
something like

.. code-block:: bash

   pileupCalc.py -i ~/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 ./2016PUHist_nominal.root

and a MC pileup profile. MC pileup profiles for official CMSSW productions are
currently hardcoded inside the ``makePUReWeightJSON.py``, and can be specified
by their tag or name in that list; the available tags can be listed by
specifying the ``--listmcprofiles`` option. The full command then becomes
something like

.. code-block:: bash

   makePUReWeightJSON.py --mcprofile "Moriond17_25ns" --nominal=2016PUHist.root --up=2016PUHist_up.root --down=2016PUHist_down.root --makePlot

To include the weight when filling plots, it is sufficient to add the weight to
a selection (usually one of the topmost in the analysis, e.g. in the
``prepareTree`` method of the analysis module).
The :py:func:`bamboo.analysisutils.makePileupWeight` method can be used to build
an expression for the weight, starting from the path of the JSON file with
weights from above, and an expression for the true number of interactions in the
event (mean of the Poissonian used), e.g. ``tree.Pileup_nTrueInt`` for NanoAOD.


.. _ugrecipetaucleaning:

Cleaning collections
''''''''''''''''''''

The CMS reconstruction sometimes ends up double-counting some objects.
This can be because of the different quality criteria used to identify each
object or because of the different performance and inner working of
the reconstruction algorithms.
Tau reconstruction for example operates on clusters that are usually
reconstructed as jets, and on top of that it can easily pick up even isolated
muons or electrons as taus (i.e. as clusters of energy with one, two, or three
tracks).

It is oftentimes necessary therefore to clean a collection of objects by
excluding any object that is spatially in the sample place of another object
whose reconstruction we trust more.

We trust more muon and electron reconstrution than tau reconstruction,
after all the quality cuts (ID efficiencies for muons and electrons are around
99.X%, whereas tau ID efficiencies are of the order of 70%.
Misidentification rates are similarly quite different), and therefore we exclude
from the tau collection any tau that happens to include within its
reconstruction cone a muon or an electron.

Bamboo provides a handy syntax for that, resulting in something like

.. code-block:: python

   cleanedTaus = op.select(taus, lambda it : op.AND(
         op.NOT(op.rng_any(electrons, lambda ie : op.deltaR(it.p4, ie.p4) < 0.3 )),
         op.NOT(op.rng_any(muons, lambda im : op.deltaR(it.p4, im.p4) < 0.3 ))
         ))

In this example, we assume that the collections ``taus``, ``electrons``, and
``muons``, have already been defined via
``taus = op.select(t.Tau, lambda tau : ...)``, and we move on to use the method
``op.rng_any()`` to filter all taus that are within a cone of a given size
(0.3, in the example) from any selected electron or muon.


.. _ugrecipejetsystematics:

Jet systematics
'''''''''''''''

For propagating uncertainties related to the jet energy calibration, and the
difference in jet energy resolution between data and simulation, each jet in
the reconstructed jet collection should be modified, the collection sorted,
and any derived quantity re-evaluated.

For efficiency and consistency, this is done by a single C++ module
that produces a set of jet collections (technically, only lists with the sorted
new momenta and the corresponding indices of the original jets are stored, and
the python decorations take care of redirecting to those if necessary).
The base jet collection proxy (e.g. ``_Jet`` for NanoAOD) has a member ``calc``
that can be used to store a reference to the module instance, and serves as a
handle to configure it.
By default, only the nominal jet collection (``Jet`` for NanoAOD) is
available. The :py:meth:`bamboo.analysisutils.configureJets` provides a
convenient way to correct the jet resolution in MC, apply a different JEC, and
add variations due to different sources of uncertainty in the jet energy scale.
As an example, the ``prepareTrees`` method of a NanoAOD analysis module could
look like this to apply a newer JEC to 2016 data and perform smearing and add
uncertainties to 2016 MC:

.. code-block:: python

   def prepareTree(self, tree, era=None, sample=None):
       ## initializes tree.Jet.calc so should be called first (better: use super() instead)
       tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, era=era, sample=sample)
       from bamboo.analysisutils import configureJets
       if era == "2016":
           if self.isMC(sample): # can be inferred from sample name
               configureJets(tree, "Jet", "AK4PFchs",
                   jec="Summer16_07Aug2017_V20_MC",
                   smear="Summer16_25nsV1_MC",
                   jesUncertaintySources=["Total"])
           else:
               if "2016G" in sample or "2016H" in sample:
                   configureJets(tree, "Jet", "AK4PFchs",
                       jec="Summer16_07Aug2017GH_V11_DATA")
               elif ...: ## other 2016 periods
                   pass

       return tree,noSel,be,lumiArgs

The jet collections ``t._Jet["nominal"]``, ``t._Jet["jerup"]``,
``t._Jet["jerdown"]``, ``t._Jet["jesTotalUp"]`` and ``t._Jet["jesTotalDown"]``
will then be available when defining plots, but they are generally not needed
directly, because by default the jet variations compatible with the
configuration will be propagated to the histograms (this can be disabled by
passing ``enableSystematics=[]`` to
:py:meth:`bamboo.analysisutils.configureJets`).

The necessary txt files will be automatically downloaded (and kept up to date)
from the repositories on github, and stored in a local cache (this should be
entirely transparent to the user |---| in case of doubt one can remove the
corresponding directory and status file from ``~/.bamboo/cache``, they will be
recreated automatically at the next use).


.. _ugreciperochester:

Rochester correction for muons
''''''''''''''''''''''''''''''

The so-called
`Rochester correction <https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon>`_
removes a bias in the muon momentum, and improves the agreement between data
and simulation in the description of the Z boson peak.
Muons are saved on the NanoAOD without this correction, so if required it needs
to be applied afterwards.
This can be done in a way that is very similar to the jet variations above,
with the difference that the module that applies the correction is attached
to the ``t._Muon`` branch, such that the corrected muons can be directly
accessed from the ``t.Muon`` branch |---| therefore the correction is
configured with

.. code-block:: python

   def prepareTree(self, tree, era=None, sample=None):
       ## initializes tree.Jet.calc so should be called first (better: use super() instead)
       tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, era=era, sample=sample)
       from bamboo.analysisutils import configureRochesterCorrection
       if era == "2016":
           configureRochesterCorrection(tree._Muon.calc, "RoccoR2016.txt")
       return tree,noSel,be,lumiArgs

The original, uncorrected, muon collection can be found in ``t._Muon.orig``.


.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _YAML: https://yaml.org

.. _YAML Wikipedia page: https://en.wikipedia.org/wiki/YAML

.. _TTree: https://root.cern/doc/master/classTTree.html

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. _DAS: https://cmsweb.cern.ch/das/

.. _SAMADhi: https://cp3.irmp.ucl.ac.be/samadhi/index.php

.. _CP3-llbb framework: https://github.com/cp3-llbb/Framework

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. _pileupcalc documentation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II

.. |---| unicode:: U+2014
   :trim:
