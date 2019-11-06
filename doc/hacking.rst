Under the hood
==============

This page collects some useful information for debugging and developing
bamboo_.

Debugging problems
------------------

Despite a number of internal checks, bamboo_ may not work correctly in some
cases.
If you encounter a problem, the following list of tips may help to find some
clues about what is going wrong:

* does the error message (python exception, or ROOT printout) provide any hint?
  (for batch jobs: check the logfile, its name should be printed)
* try to rerun on (one of) the offending sample(s), with debug printout turned
  on (by passing the ``-v`` or ``--verbose`` option, for failed batch jobs the
  main program prints the command to reproduce)
* if the problem occurs only for one or some samples: is there anything special
  in the analysis module for this sample, or in its tree format?
  The interactive mode to explore the decorated tree can be very useful to
  understand problems with expressions.
* in case of a segmentation violation while processing the events: check if you
  are not accessing any items from a container that are not guaranteed to exist
  (i.e. if you plot properties of the 2nd highest-pt jet in the event, the
  event selection should require at least two jets; with combinations or
  selections of containers this may not always be easy to find)
* check the `open issues`_ to see if your problem has already been reported, or
  is a known limitation, and, if not, ask for help on `mattermost`_ or directly
  create a `new issue`_

Different components and their interactions
-------------------------------------------

Expressions: proxies and operations
'''''''''''''''''''''''''''''''''''

The code to define expressions (cuts, weight factors, variables) in bamboo_
is designed to look as if the per-event variables (or columns of the
RDataFrame_) are manipulated directly, but what actually happens is that a
python object representation is constructed.
The classes used for this are defined in the :py:mod:`bamboo.treeoperations`
module, and inherit from :py:class:`~bamboo.treeoperations.TupleOp`.
There are currently about 25 concrete implementations.

These classes contain the minimal needed information to obtain the value they
represent (e.g. names of columns to retrieve, methods to call), but generally
no complete type information or convenience methods to use them.
They are used by almost all other bamboo_ code, but not meant to be directly
manipulated by the user code |---| this is what the proxy classes are for.

The main restriction on :py:class:`~bamboo.treeoperations.TupleOp` classes is
that, once constructed, the operation part of an expression should not be
modified.
More specifically: not after they have been passed to any backend code (so
directly after construction, e.g. by cloning, should be safe, but since
subexpressions may be passed on-demand, one should not make any assumptions in
other cases).
This allows to cache the hash of an operation, and thus very fast lookup of
expressions in sets and dictionaries, which the backend uses extensively.

The proxy classes wrap one or more operations, and behave as the resulting
value.
In some cases the correspondence is trivial, e.g. a branch with a single
floating-point number is retrieved with a
:py:class:`~bamboo.treeoperations.GetColumn` operation, and wrapped with a
:py:class:`~bamboo.treeproxies.FloatProxy`, which overloads operators for
basic math, but a proxy can also represent an object or concept that does not
correspond to a C++ type stored on the input tree, e.g. an electron (the
collection of values with the same index in all ``Electron_*[nElectron]``
branches), or a subset of the collection of electrons, whose associated
operation would be a list of indices, with the proxy holding a reference to
the original collection proxy.

All proxy classes (currently about 25) are defined in the
:py:mod:`bamboo.treeproxies` module, and inherit from the
:py:class:`~bamboo.treeoperations.TupleBaseProxy` base class, which means they
need to have an associated type, and hold a reference to a parent operation.
Operations only refer to other operations and constants, not to proxies, so
when an action (overloaded operator, member method, or a function from
:py:mod:`bamboo.treefunctions`) is performed on a proxy, a new proxy is
returned that wraps the resulting operation.

In principle proxies are only there for the user code: starting from the input
tree proxy, expressions are generated and passed to the backend, which strips
off the proxy, and generates executable code from the operation (possibly
retaining the result type from the proxy, if relevant for the produced output,
e.g. when producing a skimmed tree).
There are therefore few constraints on how the proxy classes work, as long as
the result of any action on them produces a valid operation with the expected
meaning.

Tree decorations
''''''''''''''''

All user-defined expressions start from the decorated input tree, which can,
following the previous subsection, be seen as a tree proxy.
In fact, this is exactly what the tree decoration method does: it generates the
necessary ad-hoc types that inherit from the building block proxy classes from
:py:mod:`bamboo.treeproxies`, and also have all the attributes corresponding to
the branches of the input tree.
Technically, this is done with the `type builtin`_, and a few `descriptor`_
classes.

Much of the information needed for this can be obtained by introspecting the
tree, but some details, e.g. about systematics to enable, may need to be
supplied by the user.

Selections, plots, and the RDataFrame
'''''''''''''''''''''''''''''''''''''

The main thing to know about the RDataFrame_ in bamboo_ is that partial results
are declared upon construction of :py:class:`~bamboo.plots.Plot` and
:py:class:`~bamboo.plots.Selection` objects.
The backend keeps a shadow graph of selections (with their alternatives under
systematic variations, if needed), and, for each of these, a list of the
operations that have been defined as a new column.

When an operation is converted to a C++ expression string, a reference to the
selection node where it is needed is passed, such that subexpressions can be
defined on-demand (as explained in :ref:`this section <ugcutordering>`, if a
precalculated column is needed for a selection, it may be beneficial to declare
that earlier rather than later).
This makes the verbose output a bit harder to read (to avoid redeclaring the
same function, argument names are also replaced), but ensures the correct order
of definition and reasonable efficiency.
Currently, all operations that take range arguments, and those that are
explicitly marked, are precalculated.
Function calls, notably, are not, since most are cheap to evaluate |---| this is
why expensive function calls sometimes should be explicitly requested to be
precalculated for a specific selection with
:py:meth:`bamboo.analysisutils.forceDefine`.

Organisationally, the bookkeeping code, and all the code that accesses the
interpreter and RDataFrame_ directly, is kept in
:py:mod:`bamboo.dataframebackend`, while the conversion of a
:py:class:`~bamboo.treeoperations.TupleOp` is done by its
:py:meth:`~bamboo.treeoperations.TupleOp.get_cppStr` method (many of these are
trivial, but for range-based operations, which define a helper function, they
get a bit more involved).


.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _open issues: https://gitlab.cern.ch/cp3-cms/bamboo/-/boards

.. _mattermost: https://mattermost.web.cern.ch/cms-exp/channels/bamboo

.. _new issue: https://gitlab.cern.ch/cp3-cms/bamboo/issues/new?issue%5Bassignee_id%5D=&issue%5Bmilestone_id%5D=

.. _RDataFrame: https://root.cern.ch/doc/master/classROOT_1_1RDataFrame.html

.. _type builtin: https://docs.python.org/3/library/functions.html#type

.. _descriptor: https://docs.python.org/3/reference/datamodel.html#descriptors

.. |---| unicode:: U+2014
   :trim:
