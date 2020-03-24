Recipes for common tasks
========================

.. _recipescalefactors:

Using scalefactors
------------------

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
identifier of the selection they correspond to (e.g. ``id_Loose``).
The value inside this dictionary can be either a single path to a JSON file,
or a list of ``(periods, path)`` pairs, where ``periods`` is a list of run periods, in case scalefactors for different
running periods need to be combined (the ``periods`` keyword argument to
:py:meth:`~bamboo.scalefactors.get_scalefactor` can be used to select only
a certain set of these periods).
The combination is done by either weighting or randomly sampling from the
different periods, according to the fraction of the integrated luminosity in
each (by passing ``combine="weight"`` or ``combine="sample"``, respectively).
Jet flavour tagging and dilepton (e.g. trigger) scalefactors can also be
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
     "electron_2016_94" : {
         "id_Loose"  : localize_myanalysis("Electron_EGamma_SF2D_Loose.json")
         "id_Medium" : localize_myanalysis("Electron_EGamma_SF2D_Medium.json")
         "id_Tight"  : localize_myanalysis("Electron_EGamma_SF2D_Tight.json")
     },
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
 elLooseIDSF = get_scalefactor("lepton", ("electron_2016_94", "id_Loose"), systName="elID")
 hasTwoEl = noSel.refine("hasTwoEl", cut=[ op.rng_len(electrons) > 1 ],
               weight=[ elLooseIDSF(electrons[0]), elLooseIDSF(electrons[1]) ])

 jets = op.select(t.Jet, lambda j : j.p4.Pt() > 30.)
 bJets = op.select(jets, lambda j : j.btagDeepFlavB > 0.2217) ## DeepFlavour loose b-tag working point
 deepFlavB_discriVar = { "BTagDiscri": lambda j : j.btagDeepFlavB }
 deepBLooseSF = get_scalefactor("jet", ("btag_2016_94", "DeepJet_loose"), additionalVariables=deepFlavB_discriVar, systName="bTag")
 hasTwoElTwoB = hasTwoEl.refine("hasTwoElTwoB", cut=[ op.rng_len(bJets) > 1 ],
                  weight=[ deepBLooseSF(bJets[0]), deepBLooseSF(bJets[1]) ])

Note that the user is responsible for making sure that the weights are only applied to simulated events, and not to real data!

.. _recipepureweighting:

Pileup reweighting
------------------

Pileup reweighting to make the pileup distribution in simulation match the one
in data is very similar to applying a scalefactor, except that the efficiency
correction is for the whole event or per-object |---| so the same code can be
used.
The ``makePUReWeightJSON`` script included in bamboo can be used to make
a JSON file with weights out of a data pileup profile obtained by running
``pileupcalc.py``
(inside CMSSW, see the `pileupcalc documentation`_ for details), e.g. with
something like

.. code-block:: bash

   pileupCalc.py -i ~/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt --inputLumiJSON /afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/PileUp/pileup_latest.txt --calcMode true --minBiasXsec 69200 --maxPileupBin 80 --numPileupBins 80 ./2016PUHist_nominal.root

and a MC pileup profile.
Data pileup distributions corresponding to the golden JSON files for Run 2 are
provided by the luminosity POG, see
`this hypernews annoncement<https://hypernews.cern.ch/HyperNews/CMS/get/physics-validation/3374/2.html>`_.
The MC pileup profiles for used official CMSSW productions are
currently hardcoded inside the ``makePUReWeightJSON``, and can be specified
by their tag or name in that list; the available tags can be listed by
specifying the ``--listmcprofiles`` option. The full command then becomes
something like

.. code-block:: bash

   makePUReWeightJSON --mcprofile "Moriond17_25ns" --nominal=2016PUHist.root --up=2016PUHist_up.root --down=2016PUHist_down.root --makePlot

To include the weight when filling plots, it is sufficient to add the weight to
a selection (usually one of the topmost in the analysis, e.g. in the
``prepareTree`` method of the analysis module).
The :py:func:`bamboo.analysisutils.makePileupWeight` method can be used to build
an expression for the weight, starting from the path of the JSON file with
weights from above, and an expression for the true number of interactions in the
event (mean of the Poissonian used), e.g. ``tree.Pileup_nTrueInt`` for NanoAOD.


.. _recipetaucleaning:

Cleaning collections
--------------------

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


.. _recipejetsystematics:

Jet and MET systematics
-----------------------

For propagating uncertainties related to the jet energy calibration, and the
difference in jet energy resolution between data and simulation, each jet in
the reconstructed jet collection should be modified, the collection sorted,
and any derived quantity re-evaluated.

How to do this depends on the input trees: in production NanoAOD the modified
momenta need to be calculated using the jet energy correction parameters; it is
also possible to add them when post-processing with the
`jetmetUncertainties module`_ of the NanoAODTools_ package.
In the latter case the NanoAOD decoration method will detect the modified
branches, and no additional configuration is needed.

To calculate the variations on the fly, two things are needed: when decorating
the tree, some redirections should be set up to pick up the variations from a
calculator module, and then this module needs to be configured with the correct
JEC and resolution parameters.
The first step can be done by adding ``"nJet"`` (the name of the length branch
of the stored jet collection) to the ``calcToAdd`` keyword argument of the
:py:meth:`~bamboo.NanoAODModule.prepareTree` method (which will pass this to the
:py:meth:`bamboo.treedecorators.decorateNanoAOD` method).
Similarly, to have all these variations propagated to the missing transverse
momentum, ``"MET"`` (or ``"METFixEE2017"``) should be added to the list.
Next, a calculator must be added and configured.
This can be done with the :py:meth:`bamboo.analysisutils.configureJets` and
:py:meth:`bamboo.analysisutils.configureType1MET` methods, which provide a
convenient way to correct the jet resolution in MC, apply a different JEC, and
add variations due to different sources of uncertainty in the jet energy scale,
for the jet collection and MET, respectively (the arguments should be the same
in most cases).
As an example, the relevant code of a NanoAOD analysis module could
look like this to apply a newer JEC to 2016 data and perform smearing, add
uncertainties to 2016 MC, and the same for the MET:

.. code-block:: python

   class MyAnalysisModule(NanoAODHistoModule):
       def prepareTree(self, tree, sample=None, sampleCfg=None):
           tree,noSel,be,lumiArgs = super(MyAnalysisModule, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg, calcToAdd=["nJet"])
           from bamboo.analysisutils import configureJets, configureType1MET
           isNotWorker = (self.args.distributed != "worker")
           era = sampleCfg["era"]
           if era == "2016":
               if self.isMC(sample): # can be inferred from sample name
                   configureJets(tree._Jet, "AK4PFchs",
                       jec="Summer16_07Aug2017_V20_MC",
                       smear="Summer16_25nsV1_MC",
                       jesUncertaintySources=["Total"],
                       mayWriteCache=isNotWorker,
                       isMC=self.isMC(sample), backend=be, uName=sample)
                   configureType1MET(tree._MET, "AK4PFchs",
                       jec="Summer16_07Aug2017_V20_MC",
                       smear="Summer16_25nsV1_MC",
                       jesUncertaintySources=["Total"],
                       mayWriteCache=isNotWorker,
                       isMC=self.isMC(sample), backend=be, uName=sample)
               else:
                   if "2016G" in sample or "2016H" in sample:
                       configureJets(tree._Jet, "AK4PFchs",
                           jec="Summer16_07Aug2017GH_V11_DATA",
                           mayWriteCache=isNotWorker,
                           isMC=self.isMC(sample), backend=be, uName=sample)
                       configureType1MET(tree._MET, "AK4PFchs",
                           jec="Summer16_07Aug2017GH_V11_DATA",
                           mayWriteCache=isNotWorker,
                           isMC=self.isMC(sample), backend=be, uName=sample)
                   elif ...: ## other 2016 periods
                       pass

           return tree,noSel,be,lumiArgs

Both with variations read from a postprocessed NanoAOD and calculated on the
fly, the different jet collections are available from ``t._Jet``, e.g.
``t._Jet["nom"]`` (postprocessed) or ``t._Jet["nominal"]`` (calculated),
``t._Jet["jerup"]``, ``t._Jet["jerdown"]``, ``t._Jet["jesTotalUp"]``,
``t._Jet["jesTotalDown"]`` etc. depending on the configured variations
(when accessing these directly, ``t._Jet[variation][j.idx]`` should be used
to retrieve the entry corresponding to a specific jet ``j``, if the latter is
obtained from a selected and/or sorted version of the original collection |---|
``object.idx`` is always the index in the collection as found in the tree).

``t.Jet`` will be changed for one of the above for each systematic variation,
if it affects a plot, in case automatically producing the systematic variations
is enabled (the collections from ``t._Jet`` will not be changed).
The automatic calculation of systematic variations can be disabled globally
or on a per-selection basis (see above), and for on the fly calculation also by
passing ``enableSystematics=[]`` to
:py:meth:`bamboo.analysisutils.configureJets`).
The jet collection as stored on the input file, finally, can be retrieved as
``t._Jet.orig``.

.. important:: Sorting the jets
   No sorting is done as part of the above procedure, so if relevant this
   should be added by the user (e.g. using
   ``jets = op.sort(t.Jet, lambda j : -j.pt)`` for sorting by decreasing
   transverse momentum).
   In a previous version of the code this was included, but since some selection
   is usually applied on the jets anyway, it is simpler (and more efficient) to
   perform the sorting then.

.. note:: Isn't it slow to calculate jet corrections on the fly?
   It does take a bit of time, but the calculation is done in one C++ module,
   which should not be executed more than once per event (see the explanation
   of the :py:meth:`bamboo.analysisutils.forceDefine` method in the
   :ref:`section above<ugcutordering>`).
   In most realistic cases, the bottleneck is in reading and decompressing the
   input files, so the performance hit from the jet corrections should usually
   be acceptable.

.. tip:: Bamboo_ runs outside CMSSW and has no access to the conditions
   database, so it fetches the necessary txt files from the repositories
   on github (they are quite large, so this is more efficient than storing
   a clone locally). They should automatically be updated if the upstream
   repository changes and the ``mayWriteCache`` argument to
   :py:meth:`bamboo.analysisutils.configureJets` (see the example above)
   helps ensure that only one process write to the cache at a time.
   In case of doubt one can use the ``checkBambooCMSJetDatabaseCaches`` script
   to update or check the cache interactively and, as a last resort, remove
   the cache directories and status files from ``~/.bamboo/cache``:
   they will be recreated automatically at the next use.


.. _reciperochester:

Rochester correction for muons
------------------------------

The so-called
`Rochester correction <https://twiki.cern.ch/twiki/bin/viewauth/CMS/RochcorMuon>`_
removes a bias in the muon momentum, and improves the agreement between data
and simulation in the description of the Z boson peak.
As for the jet correction and variations described in the previous section,
this can either be done during postprocessing, with the
`muonScaleResProducer module`_ of the NanoAODTools_ package, or on the fly.
The former will be detected automatically when decorating the trees, the latter
requires adding ``"nMuon"`` to the ``calcToAdd`` keyword argument passed to he
:py:attr:`~bamboo.NanoAODModule.prepareTree` method, such that the decorations
will pick up the momenta from a calculator module
instead of directly from the input file.

The on the fly calculator should be added and configured with the
:py:meth:`bamboo.analysisutils.configureRochesterCorrection` method,
as in the example below.
``tree._Muon`` keeps track of everything related to the calculator; the
uncorrected muon collection can be found in ``tree._Muon.orig``, and the
corrected muons are in ``tree.Muon``.

.. code-block:: python

   class MyAnalysisModule(NanoAODHistoModule):
       def prepareTree(self, tree, sample=None, sampleCfg=None):
           tree,noSel,be,lumiArgs = NanoAODHistoModule.prepareTree(self, tree, sample=sample, sampleCfg=sampleCfg, calcToAdd=["nMuon"])
           from bamboo.analysisutils import configureRochesterCorrection
           era = sampleCfg["era"]
           if era == "2016":
               configureRochesterCorrection(tree._Muon, "RoccoR2016.txt", isMC=self.isMC(sample), backend=be, uName=sample)
       return tree,noSel,be,lumiArgs

.. _recipesplitsamplesubcomp:

Splitting a sample into sub-components
--------------------------------------

It is frequently necessary to split a single Monte-Carlo sample into different processes, depending on generator-level information, or simply to add some cuts at generator level (e.g. to stitch binned samples together). 
This can be achieved by duplicating that sample in the analysis configuration file for as many splits as are needed, and putting (any) additional information into that sample's entry, e.g. as:

.. code-block:: yaml

     ttbb:
       db: das:/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM
       era: 2018
       group: ttbb
       subprocess: ttbb
       cross-section: 366.
       generated-events: genEventSumw

     ttjj:
       db: das:/TTToSemiLeptonic_TuneCP5_13TeV-powheg-pythia8/RunIIAutumn18NanoAODv5-Nano1June2019_102X_upgrade2018_realistic_v19-v1/NANOAODSIM
       era: 2018
       group: ttjj
       subprocess: ttjj
       cross-section: 366.
       generated-events: genEventSumw

That information can then be retrieved in the analysis module through the ``sampleCfg`` keyword argument, to add additional cuts to the selection when preparing the tree:

.. code-block:: python

   def prepareTree(self, tree, sample=None, sampleCfg=None):
       tree,noSel,be,lumiArgs = super(MyAnalysisModule, self).prepareTree(tree, sample=sample, sampleCfg=sampleCfg)

       if "subprocess" in sampleCfg:
            subProc = sampleCfg["subprocess"]
            if subProc == "ttbb":
                noSel = noSel.refine(subProc, cut=(tree.genTtbarId % 100) >= 52)
            elif subProc == "ttjj":
                noSel = noSel.refine(subProc, cut=(tree.genTtbarId % 100) < 41)

       return tree,noSel,be,lumiArgs


.. _recipecmdlinearg:

Adding command-line arguments
-----------------------------

The base :ref:`analysis module<uganalysismodule>`,
:py:class:`bamboo.analysismodules.AnalysisModule`, calls the
:py:meth:`~bamboo.analysismodules.AnalysisModule.addArgs` method (the default
implementation does nothing) when constructing the command-line arguments
parser (using the `argparse`_ module).
Analysis modules can reimplement this method to specify more arguments, e.g.

.. code-block:: python

    class MyModule(...):

        def addArgs(self, parser):
            super(MyModule, self).addArgs(parser)
            parser.add_argument("--whichPlots", type=str,
                                default="control",
                                help="Set of plots to produce")


The parsed arguments are available under the ``args`` member variable, e.g.
``self.args.whichPlots`` in the example above.
The complete list of command-line options (including those specified in the
analysis module) can be printed with ``bambooRun -h -m myModule.py.MyModule``.
In fact the parser argument is an
`argument group`_,
so they are listed separately from those in the base class.
This is also used to copy all user-defined arguments to the commands that are
passed to the worker tasks, when running in distributed mode.

.. _recipecustomanacfg:

Editing the analysis configuration
----------------------------------

Similarly to the above, it is possible to modify the analysis configuration 
(loaded from the YAML file) from a module before the configuration 
is used to create jobs (in distributed mode), run on any file (in sequential mode),
or run plotIt (in the postprocessing step).
This allows e.g. to change the samples that are going to be used, change the list 
of systematics, etc., without having to edit manually the YAML file or maintaining separate files.
Below is an example of how this works:

.. code-block:: python

    class MyModule(...):

        def customizeAnalysisCfg(self, analysisCfg):
            for smp in list(analysisCfg["samples"]):
                if not analysisCfg["samples"][smp].get("is_signal", False):
                    del analysisCfg["samples"][smp]
            


.. _recipemvaevaluate:

Evaluate an MVA classifier
--------------------------

This feature is still under development, a recipe will be added here shortly

.. _recipemergedcategoryplots:

Make combined plots for different selections
--------------------------------------------

It is rather common to define categories with e.g. different lepton flavours
and selections, but then make plots with the entries from these (disjoint)
sets of events combined.
Given the structure of the RDataFrame_ graph and the
:py:class:`~bamboo.plots.Selection` tree, the most convenient way to achieve
this is by defining the histograms for each category, and make a merged
histogram later on.
The :py:class:`~bamboo.plots.SummedPlot` class does exactly this, and since it
presents the same interface to the analysis module as a regular
:py:class:`~bamboo.plots.Plot`, it can simply be added to the same plot list
(to produce only the combined plot and not those for the individual
contributions, it is sufficient to not add the latter to the plot list), e.g.

.. code-block:: python

   from bamboo.plots import Plot, SummedPlot, EquidistantBinning
   mjj_mumu = Plot.make1D("Mjj_MuMu", op.invariant_mass(jets[0].p4, jets[1].p4),
                          sel_mumu, EquidistantBinning(50, 20., 120.))
   mjj_elel = Plot.make1D("Mjj_ElEl", op.invariant_mass(jets[0].p4, jets[1].p4),
                          sel_elel, EquidistantBinning(50, 20., 120.))
   mjj_sum = SummedPlot("Mjj", [mjj_mumu, mjj_elel], title="m(jj)")
   plots += [ mjj_mumu, mjj_elel, mjj_sum ] # produce all plots


The other plot properties of a :py:class:`~bamboo.plots.SummedPlot` (titles,
labels etc.) can be specified with keyword arguments to the constructor;
otherwise they are taken from the first component plot.

.. note:: :py:class:`~bamboo.plots.SummedPlot` simply adds up the histograms,
   it is up to the user to make sure an event can only enter one of the
   categories, if this is what it is used for.


.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _CP3-llbb framework: https://github.com/cp3-llbb/Framework

.. _pileupcalc documentation: https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Pileup_JSON_Files_For_Run_II

.. _NanoAODTools: https://github.com/cms-nanoAOD/nanoAOD-tools

.. _jetmetUncertainties module: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/jme/jetmetUncertainties.py

.. _muonScaleResProducer module: https://github.com/cms-nanoAOD/nanoAOD-tools/blob/master/python/postprocessing/modules/common/muonScaleResProducer.py

.. _argparse: https://docs.python.org/3/library/argparse.html

.. _argument group: https://docs.python.org/3/library/argparse.html#argument-groups

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. |---| unicode:: U+2014
   :trim:
