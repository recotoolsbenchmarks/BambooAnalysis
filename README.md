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

## Documentation

The HTML documentation (with a longer introduction, installation instructions,
recipes for common tasks and an API reference of the classes and methdos) is
available [here](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/)
or can be generated from the source code repository with
```sh
python setup.py build
python setup.py build_sphinx
```
(the output will be placed in `doc/build/html`).

## Development

Bamboo is currently at the beta stage, which means that a first version of most
of the planned features is there and it starts to be usable, but some things are
missing, incomplete, or buggy - your feedback is essential for fixing those, so
please watch out, check, feel free to add tests if something is not covered yet,
and do not hesitate to
[open an issue](https://gitlab.cern.ch/cp3-cms/bamboo/boards) or ask in
[the mattermost channel](https://cp3-mm.irmp.ucl.ac.be/cp3-llbb/channels/bamboo)
in case of doubt.
