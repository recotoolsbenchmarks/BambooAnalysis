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


Tree decorations
''''''''''''''''


Selections, plots, and the RDataFrame
'''''''''''''''''''''''''''''''''''''



.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _open issues: https://gitlab.cern.ch/cp3-cms/bamboo/-/boards

.. _mattermost: https://mattermost.web.cern.ch/cms-exp/channels/bamboo

.. _new issue: https://gitlab.cern.ch/cp3-cms/bamboo/issues/new?issue%5Bassignee_id%5D=&issue%5Bmilestone_id%5D=
