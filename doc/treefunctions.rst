Building expressions
====================

In order to efficiently process the input files, bamboo builds up an object representation
of the expressions (cuts, weights, plotted variables) needed to fill the histograms, and
dynamically generates C++ code that is passed to RDataFrame_.
The expression trees are built up throug proxy classes, which mimic the final type (there are
e.g. integer and floating-point number proxy classes that overload the basic mathematical operators),
and generate a new proxy when called.
As an example: ``t.Muon[0].charge`` gives an integer proxy to the operation corresponding
to ``Muon_charge[0]``; when the addition operator is called in ``t.Muon[0].charge+t.Muon[1].charge``,
an integer proxy to (the object representation of) ``Muon_charge[0]+Muon_charge[1]`` is returned.

The proxy classes try to behave as much as possible as the objects they represent, so in most cases
they can be used as if they really were a number, boolean, momentum fourvector... or a muon,
electron, jet etc. |---| simple 'struct' types for those are generated when decorating the tree,
based on the branches that are found.
Some operations, however, cannot easily be implemented in this way, for instance mathematical functions and
operations on containers. Therefore, the :py:mod:`bamboo.treefunctions` module provides a set of
additional helper methods ( such that the user does not need to know about the implementation details
in the :py:mod:`bamboo.treeoperations` and :py:mod:`bamboo.treeproxies` modules). In order to keep
the analysis code compact, it is recommended to import it with

.. code-block:: python

   from bamboo import treefunctions as op 

inside every analysis module. The available functions are listed below.

List of functions
-----------------

.. automodule:: bamboo.treefunctions
   :members:
   :member-order: bysource

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. |---| unicode:: U+2014
   :trim:
