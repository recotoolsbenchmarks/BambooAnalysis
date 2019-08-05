Bamboo: A high-level HEP analysis library for ROOT::RDataFrame
==============================================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   install
   userguide
   treefunctions
   apiref

The RDataFrame_ class provides an efficient and flexible way to process
per-event information (stored in a TTree_) and e.g. aggregate it into
histograms.

With the typical pattern of storing object arrays as a structure of arrays
(variable-sized branches with a common prefix in the names and length),
the expressions that are typically needed for a complete analysis quickly become
cumbersome to write (with indices to match, repeated sub-expressions etc.).
As an example, imagine the expression needed to calculate the invariant mass of
the two leading muons from a NanoAOD (which stores momenta with pt, eta and phi
branches |---| one way is to construct LorentzVector objects, sum and evaluate
the invariant mass.
Next imagine doing the same thing with the two highest-pt jets that have a b-tag
and are not within some cone of the two leptons you already selected in another
way (while keeping the code maintainable enough to allow for passing jet momenta
with a systematic variation applied).

Bamboo attempts to solve this problem by automatically constructing
lightweight python wrappers based on the structure of the TTree_,
which allow to construct such expression with high-level code, similar to the
language that is commonly used to discuss and describe them. By constructing
an object representation of the expression, a few powerful operations can be
used to compose complex expressions.
This also allows to automate the construction of derived expressions, e.g. for
shape systematic variation histograms.

Building selections, plots etc. with such expressions is analysis-specific, but
the mechanics of loading data samples, processing them locally or on a batch
system (and merging the output of that), combining the outputs for different
samples in an overview etc. is very similar over a broad range of use cases.
Therefore a common implementation of these is provided, which can be used by
extending a base class (to fill histograms and make stacked plots, a class
needs to be written with a method that returns a list of 'plot' objects |---|
each essentially a combination of an x-axis variable, selection, and weight
to apply to every event |---| and a configuration file that specifies which
datasets should be processed, and how they should be stacked).

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. _TTree: https://root.cern/doc/master/classTTree.html

.. |---| unicode:: U+2014
   :trim:

Indices and tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
