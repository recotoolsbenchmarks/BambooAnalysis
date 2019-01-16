# Bamboo: A high-level HEP analysis library for ROOT::RDataFrame

The [`ROOT::RDataFrame`](https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html)
class provides an efficient and flexible way to process per-event information
(stored in a `TTree`) and e.g. aggregate it into histograms.
With the typical pattern of storing object arrays as a structure of arrays
(variable-sized branches with a common prefix in the names and length),
the expressions that are typically needed for a complete analysis quickly become
cumbersome to write, with indices to match, repeated sub-expressions etc. -

Bamboo attempts to solve this problem by automatically constructing
lightweight python wrappers based on the structure of the `TTree`,
which allow to construct such expression with high-level code, similar to the
language that is commonly used to discuss and describe them. By constructing
an object representation of the expression, a few powerful operations can be
used to compose complex expressions.

Building selections, plots etc. with such expressions is analysis-specific, but
the mechanics of loading data samples, processing them locally or on a batch
system, combining the outputs for different samples in an overview etc.
is very similar over a broad range of use cases.
Therefore a common implementation of these is provided, such that the analyst
only needs to provide a subclass with their selection and plot definitions,
and a configuration file with a list of samples, and instructions how to
display them.

To view the documentation, for now (i.e. until automatic deployment to a website
is enabled) run
```sh
python setup.py build ## only needed if you haven't installed bamboo yet
python setup.py build_sphinx
```
inside your bamboo clone. This will generate the documentation in
`doc/build/html`.
