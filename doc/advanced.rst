Advanced topics
===============

.. _ugcppmodules:

Loading (and using) C++ modules
-------------------------------

The :py:mod:`bamboo.root` module defines a few thin wrappers to load additional
libraries or headers with the ROOT interpreter (PyROOT makes them directly
accessible through the global namespace, which can be imported as ``gbl`` from
:py:mod:`bamboo.root`.
As an example, the library that contains the dictionaries for the classes
used in Delphes_ output trees can be added with the
:py:meth:`bamboo.root.loadLibrary` method:

.. code-block:: python

    from bamboo.root import loadLibrary
    loadLibrary("libDelphes")

In this specific case, this line should be added to the analysis module's
:py:meth:`~bamboo.analysismodules.AnalysisModule.initialize` method, to make
sure the library is added before any file is opened.
For a module that is used in calculating expressions, it is sufficient to load
it from the :py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree` or
:py:meth:`~bamboo.analysismodules.HistogramsModule.processTrees` method
(repeated loads should not cause problems).

C++ methods can be used directly from an expression.
Non-member methods that are known by the interpreter (e.g. because the
corresponding header has been added with :py:meth:`bamboo.root.loadHeader`),
can be retrieved with the :py:meth:`bamboo.treefunctions.extMethod`.

It is often useful to define a class that stores some parameters, and then call
a member method with event quantities to obtain a derived quantity (this is
also the mechanism used for most of the builtin corrections).
In order to use such a class, its header (and shared library, if necessary)
should be loaded as above, and an instance defined with
:py:meth:`bamboo.treefunctions.define`, e.g.

.. code-block:: python

    myCalc = op.define("MyCalculatorType", 'const auto <<name>> = MyCalculatorType("test");')
    myCorrection = myCalc.evaluate(tree.Muon[0].pt, tree.Muon[1].pt)

.. warning:: With implicit multi-threading enabled, only thread-safe methods may
    be called in this way (e.g. const member methods, without global or member
    variables used for caching).

.. note:: The usual logic to avoid redefinition of these variables is applied.
    In cases like above where all parameters are supplied at once, this will
    work as expected.
    If the calculator is further configured by calling member methods (it can
    be accessed directly through PyROOT), it is safer to create a unique
    instance for each sample, e.g. by adding a comment that contains the sample
    name at the end of the declaration (an optional ``nameHint`` argument can
    be given to make the generated code more readable, but this will be ignored
    in case the declaration string is the same).

.. _ugcutordering:

Ordering selections and plots efficiently
-----------------------------------------

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

   def definePlots(self, t, noSel, sample=None, sampleCfg=None):
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


.. _Delphes: https://cp3.irmp.ucl.ac.be/projects/delphes

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. |---| unicode:: U+2014
   :trim:
