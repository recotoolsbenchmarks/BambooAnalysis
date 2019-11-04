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
* the ``--eras`` option specifies which of the eras from the configuration file
  to consider, and which type of plots to make. The format is
  ``[mode][:][era1,era2,...]``, where ``mode`` is one of ``split`` (plots for
  each of the eras separately), ``combined`` (only plots for all eras combined)
  or ``all`` (both of these, this is the default).

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

Three top-level keys are currently required: ``tree`` with the name of the TTree_
inside the file (e.g. ``tree: Events`` for NanoAOD), ``samples`` with a list of
samples to consider, and ``eras``, with a list of data-taking periods and their
integrated luminosity.
For stacked histogram plots, a ``plotIt`` section should also be specified (the
:py:func:`bamboo.analysisutils.runPlotIt` method will insert the ``files`` and
``plots`` sections and run plotIt_ with the resulting configuration; depending
on the ``--eras`` option passed, per-era or combined plots will be produced, or
both, which is the default).

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
The ``era`` key specifies which era (one of those specified in the ``eras``
section, see above) the sample corresponds to, and which luminosity value
should be used for the normalisation.

For the normalization of simulated samples in the stacks, the number of
generated evens and cross-section are also needed. The latter should be
specified as ``cross-section`` with the sample (in the same units as the
``luminosity`` for the corresponding ``era``), the former can be computed from
the input files. For this, the
:py:class:`bamboo.analysismodules.HistogramsModule` base class will call the
``mergeCounters`` method when processing the samples, and the ``readCounters``
method to read the values from the results file - for NanoAOD the former merges
the `Runs` trees and saves the results, while the latter performs the sum of
the branch with the name specified under ``generated-events``.

For large samples, a ``split`` property can be specified, such that the input
files are spread out over different batch jobs.
A positive number is taken as the number of jobs to divide the inputs over,
while a negative number gives the number of files per job (e.g. ``split: 3``
An ``era`` key is also foreseen (to make 2016/2017/2018/combined plots) - but
it is currently ignored.
will create three jobs to process the sample, while ``split: -3`` will result
in jobs that process three files each).

All together a typical analysis YAML_ file would look like the following (but
with many more sample blocks, and typically a few era blocks; the ``plotIt``
section is left out for brevity).

.. code-block:: yaml

     tree: Events
     eras:
       '2016':
         luminosity: 12.34
     samples:
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

       def definePlots(self, t, noSel, sample=None, sampleCfg=None):
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



.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _YAML: https://yaml.org

.. _YAML Wikipedia page: https://en.wikipedia.org/wiki/YAML

.. _TTree: https://root.cern/doc/master/classTTree.html

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. _DAS: https://cmsweb.cern.ch/das/

.. _SAMADhi: https://cp3.irmp.ucl.ac.be/samadhi/index.php

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. |---| unicode:: U+2014
   :trim:

