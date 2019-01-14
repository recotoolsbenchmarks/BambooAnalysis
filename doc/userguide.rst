User guide
==========

This section contains some more information on doing your analysis with bamboo.
It assumes you have successfully installed it following the instructions in the
previous section.

The first thing to make sure is that bamboo can work with your trees.
A first set of decorators for CMS NanoAOD is included in
:py:func:`bamboo.treedecorators.decorateNanoAOD`; to make stacked histogram
plots from them it is sufficient to make your analysis module inherit from
:py:class:`bamboo.analysismodules.NanoAODHistoModule`
(which calls this method from its
:py:meth:`~bamboo.analysismodules.NanoAODHistoModule.prepareTrees` method).
Other types of trees can be included in a similar way, but a bit of development
is needed to provided a convenient way to do so (help welcome).

Then, you will need to provide the two really analysis-specific parts: a
configuration file that lists the samples and plot configuration, and a module
specifying what to fill in the histograms. Let's start with the latter.

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
methods that should be further specified by subclasses.

:py:class:`~bamboo.analsysimodules.HistogramsModule` does this for the
stacked histogram plots, composing
:py:meth:`~bamboo.analsysimodules.HistogramsModule.processTrees` from
:py:meth:`~bamboo.analsysimodules.HistogramsModule.prepareTree` and
:py:meth:`~bamboo.analsysimodules.HistogramsModule.definePlots`, while taking
the JSON lumi block mask and counter merging into account.
It also calls the `plotIt` executable from 
:py:meth:`~bamboo.analsysimodules.HistogramsModule.postProcess` (with the plots
list and analysis configuration file, it has all required information for that).
:py:class:`~bamboo.analysismodules.NanoAODHistoModule` supplements this with
the decorations and counter merging and reading for NanoAOD,
such that all the final module needs to do is defining plots and selections,
as in the example :py:mod:`examples.nanozmumu`.
This layered structure is used such that code can be maximally reused for other
types of trees.

For the code inside the module, the example is also very instructive:

.. code-block:: python

       def definePlots(self, t, noSel, systVar="nominal"):
           from bamboo.plots import Plot, EquidistantBinning
           from bamboo import treefunctions as op

           plots = []

           twoMuSel = noSel.refine("twoMuons", cut=[ op.rng_len(t.Muon) > 1 ])
           plots.append(Plot.make1D("dimu_M", op.invariant_mass(t.Muon[0].p4, t.Muon[1].p4), twoMuSel,
                   EquidistantBinning(100, 20., 120.), title="Dimuon invariant mass", plotopts={"show-overflow":False}))

           return plots

The key classes are defined in :py:mod:`bamboo.plots`:
:py:class:`~bamboo.plots.Plot` and :py:class:`~bamboo.plots.Selection`.
The latter represents a consistent set of selection requirements (cuts) and
weight factors (e.g. to apply corrections). Selections are defined by refining
a "root selection" with additional cuts and weights, and each should have a
unique name (an exception is raised at construction otherwise).
The root selection allows to do some customisation upfront, e.g. the applying
the JSON luminosity block mask for data. A plot object refers to a selection,
and specifies which variable(s) to plot, with which binning(s), labels, options
etc. (the ``plotOpts`` dictionary is copied directly into the plot section of the
plotIt configuration file).

Specifying cuts, weight, and variables: expressions
---------------------------------------------------

The first argument to the
:py:meth:`~bamboo.analsysimodules.HistogramsModule.definePlots`
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
:py:mod:`bamboo.treefunctions` module, so it is recommended to

.. code-block:: python

   from bamboo import treefunctions as op

(or any other name, as you prefer) in your analysis module.

Ideally, the decorated tree and the :py:mod:`bamboo.treefunctions` module
are all you ever need to import and know about the decorations.
Therefore the best way to proceed now is get a decorated tree
inside an IPython shell and play around.
For :py:mod:`bamboo.analysismodules.HistogramsModule` this can always be done
by passing the ``--interactive`` flag, with either one of
(depending on if you copied the NanoAOD test file above)

.. code-block:: sh

   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive --distributed=worker /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root
   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu --interactive bamboo/examples/test_nanozmm.yml [ --envConfig=bamboo/examples/ingrid.ini ] -o int1

The decorated tree is in the ``tree`` variable (the original ``TChain`` is in
``tup``) and the :py:mod:`bamboo.treefunctions` module is there as `op`
(the ``c_...`` methods construct a constant, whereas the ``rng_...`` methods
work on a collection and return a single value
(note: all these methods should be documented),
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

.. _YAML: https://yaml.org

.. _YAML Wikipedia page: https://en.wikipedia.org/wiki/YAML

.. _TTree: https://root.cern/doc/master/classTTree.html

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. _DAS: https://cmsweb.cern.ch/das/

.. _SAMADhi: https://cp3.irmp.ucl.ac.be/samadhi/index.php

.. |---| unicode:: U+2014
   :trim:
