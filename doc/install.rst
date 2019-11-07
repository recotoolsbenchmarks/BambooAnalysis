Installation and setup
======================

Dependencies and environment
----------------------------

Bamboo_ only depends on python3 (with pip/setuptools to install PyYAML and
numpy if needed) and a recent version of ROOT (6.16.00 is recommended because
it brings some changes to RDataFrame with respect to 6.14.06).
On ingrid and lxplus (or any machine with cvmfs), an easy way to get such
a recent version of ROOT is through a CMSSW release that depends on it (``10_4``
has a recent ROOT, but no python3 yet), or from the lcgsoft distribution, e.g.

.. code-block:: sh

   source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
   python -m venv bamboovenv
   source bamboovenv/bin/activate

(the second command creates a `virtual environment`_
to install python packages in, after installation it is sufficient to run two
other commands, to pick up the correct base system and then the installed
packages).

Some features bring in additional dependencies. Bamboo_ should detect if these
are relied on and missing, and print a clear error message.
Currently, they include:

- the dasgoclient executable (and a valid grid proxy) for retrieving the list
  of files in samples specified with ``db: das:/X/Y/Z``. Due to some
  interference with the setup script above, the best is to run the cms
  environment scripts first, and also run ``voms-proxy-init`` then (this can
  alternatively also be done from a different shell on the same machine)
- the SAMADhi_ python library, to retrieve the list of files in samples
  specified with ``db: SAMADhi:<ID-or-name>`` (not yet integrated)
- the slurm command-line tools, and CP3SlurmUtils_, which can be loaded with
  ``module load slurm/slurm_utils`` on the ingrid ui machines
- the ``makePUReWeightJSON`` script can make a plot of the PU distributions
  and weights if the ``--makePlot`` option is given, but it needs matplotlib_
  for that

Installation
------------

Bamboo_ can (and should, in most cases) be installed in a
`virtual environment`_ (see above) with pip, minimally either one of

.. code-block:: sh

   pip install git+https://gitlab.cern.ch/cp3-cms/bamboo.git
   pip install git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git

but in the current development stage it may be useful to install from
a local clone, such that you can use it to test and propose changes, using

.. code-block:: sh

   git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git /path/to/your/bambooclone
   pip install /path/to/your/bambooclone ## e.g. ./bamboo (not bamboo - a package with that name exists)

such that you can update later on with (inside ``/path/to/your/bambooclone``)

.. code-block:: sh

   git pull upstream master
   pip install --upgrade .

(installing in editable mode does not work because some C++ libraries that are
built as extensions are not correctly picked up then).

When developing bamboo, setuptools can be used to build this documentation and
run the tests (the final clean command is to prevent the next upgrade with pip
from picking up the wrong build directory).

.. code-block:: sh

   python setup.py build
   python setup.py build_sphinx
   python setup.py test
   python setup.py clean --all ## before upgrading the virtualenv with pip

.. note::

   bamboo is a shared package, so everything that is specific to a single
   analysis (or a few related analyses) is best stored elsewhere (e.g. in
   ``bamboodev/myanalysis`` in the example below); otherwise you will need to
   be very careful when updating to a newer version.

   The ``bambooRun`` command can pick up code in different ways, so it is
   possible to start from a single python file, and move to a pip-installed
   analysis package later on when code needs to be shared between modules.

For combining the different histograms in stacks and producing pdf or png files,
which is used in many analyses, the plotIt_ tool is used.
It can be installed with cmake, e.g.

.. code-block:: sh

   git clone -o upstream https://github.com/cp3-llbb/plotIt.git /path/to/your/plotitclone
   mkdir build-plotit
   cd build-plotit
   cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV /path/to/your/plotitclone
   make -j2 install
   cd -

where ``-DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV`` ensures that ``make install``
will put the ``plotIt`` executable directly in the ``bin`` directory of the
virtualenv (if not using a virtualenv, its path can be passed to ``bambooRun``
with the ``--plotIt`` command-line option).


For the impatient: recipes for installing and updating
''''''''''''''''''''''''''''''''''''''''''''''''''''''

Putting the above commands together, the following should give you a virtual
environment with bamboo_, and a clone of bamboo_ and plotIt in case you need to
modify them, all under ``bamboodev``:

Fresh install
#############

.. code-block:: sh

   mkdir bamboodev
   cd bamboodev
   # make a virtualenv
   source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
   python -m venv bamboovenv
   source bamboovenv/bin/activate
   # clone and install bamboo
   git clone -o upstream https://gitlab.cern.ch/cp3-cms/bamboo.git
   pip install ./bamboo
   # clone and install plotIt
   git clone -o upstream https://github.com/cp3-llbb/plotIt.git
   mkdir build-plotit
   cd build-plotit
   cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
   make -j2 install
   cd -

Environment setup
#################

Once bamboo_ and plotIt have been installed as above, only the following two
commands are needed to set up the environment in a new shell:

.. code-block:: sh

   source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
   source bamboodev/bamboovenv/bin/activate

Update bamboo
#############

Assuming the environment is set up as above; this can also be used to test a
pull request or local modifications to the bamboo_ source code

.. code-block:: sh

   cd bamboodev/bamboo
   git checkout master
   git pull upstream master
   pip install --upgrade .

Update plotIt
#############

Assuming the environment is set up as above; this can also be used to test a
pull request or local modifications to the plotIt source code.
If a plotIt build directory already exists it should have been created with the same
environment, otherwise the safest solution is to remove it.

.. code-block:: sh

   cd bamboodev
   mkdir build-plotIt
   cd build-plotit
   cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
   make -j2 install
   cd -

Move to a new LCG release or install an independent version
############################################################

Different virtual environments can exist alongside each other, as long as for
each the corresponding base LCG distribution is setup in a fresh shell.
This allows to have e.g. one stable version used for analysis, and another one
to test experimental changes, or check a new LCG release, without touching a
known working version.

.. code-block:: sh

   cd bamboodev
   source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
   python -m venv bamboovenv_X
   source bamboovenv_X/bin/activate
   pip install ./bamboo
   # install plotIt (as in "Update plotIt" above)
   mkdir build-plotit
   cd build-plotit
   cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
   make -j2 install
   cd -

Test your setup
---------------

Now you can run a few simple tests on a CMS NanoAOD to see if the installation
was successful. A minimal example is run by the following command:

.. code-block:: sh

   bambooRun -m /path/to/your/bambooclone/examples/nanozmumu.py:NanoZMuMu /path/to/your/bambooclone/examples/test1.yml -o test1

which will run over a single sample of ten events and fill some histograms
(in fact, only one event passes the selection, so they will not look very
interesting).
If you have a NanoAOD file with muon triggers around, you can put its path
instead of the test file in the yml file and rerun to get a nicer plot (xrootd
also works, but only for this kind of tests |---| in any practical case the
performance benefit of having the files locally is worth the cost of replicating
them).

Getting started
---------------

The test command above shows how bamboo is typically run: using the
:ref:`bambooRun<ugbambooRun>` command, with a python module that specifies what
to run, and an :ref:`analysis YAML file<uganalysisyaml>` that specifies which
samples to process, and how to combine them in plots (there are several options
to run a small test, or submit jobs to the batch system when processing a lot
of samples).

A more realistic analysis YAML configuration file is
`bamboo/examples/analysis_zmm.yml <https://gitlab.cern.ch/cp3-cms/bamboo/blob/master/examples/analysis_zmm.yml>`_,
which runs on a significant fraction of the 2016 and 2017 ``DoubleMuon`` data
and the corresponding Drell-Yan simulated samples.
Since the samples are specified by their DAS path in this case, the
``dasgoclient`` executable and a valid grid proxy are needed for resolving
those to files, and a :ref:`configuration file<ugenvconfig>` that describes the
local computing environment (i.e. the root path of the local CMS grid storage,
or the name of the redirector in case of using xrootd); examples are included
for the UCLouvain-CP3 and CERN (lxplus/lxbatch) cases.

The corresponding
`python module <https://gitlab.cern.ch/cp3-cms/bamboo/blob/master/examples/nanozmumu.py>`_
shows the typical structure of ever tighter event selections that derive from
the base selection, which accepts all the events in the input, and plots that
are defined based on these selection, and returned in a list from the main
method (this corresponds to the pdf or png files that will be produced).

The module deals with a decorated version of the tree, which can also be
inspected from an IPython shell by using the ``-i`` option to ``bambooRun``,
e.g.

.. code-block:: sh

   bambooRun -i -m /path/to/your/bambooclone/examples/nanozmumu.py:NanoZMuMu /path/to/your/bambooclone/examples/test1.yml

together with the helper methods defined on :doc:`this page<treefunctions>`,
this allows to define a wide variety of selection requirements and variables.

The :doc:`user guide<userguide>` contains a much more detailed description of
the different files and how they are used, and the
:doc:`analysis recipes page<recipes>` provides more information about the
bundled helper methods for common tasks.
The :doc:`API reference<apiref>` describes all available user-facing methods
and classes.
If the builtin functionality is not sufficient, some hints on extending or
modifying bamboo can be found in the :doc:`advanced topics<advanced>` and the
:doc:`hacking guide<hacking>`.

Machine learning packages
-------------------------

In order to evaluate machine learning classifiers, bamboo_ needs to find the
necessary C(++) libraries, both when the extension libraries are compiled and
at runtime (so they need to be installed before (re)installing bamboo_).
For libtorch_ this is done by searching for the ``torch`` package
using ``pkg_resources``, which should work whenever it is installed with pip.
For Tensorflow-C_ there is currently no pip-installable package, so it will be
searched for (by cmake and the dynamic library loader) in the default locations,
supplemented with the currently active `virtual environment`_ (a script to
install it directly there is bundled with the bamboo_ source distribution, in
``ext/install_tensorflow-c.sh``).

.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _CP3SlurmUtils: https://cp3-git.irmp.ucl.ac.be/cp3-support/helpdesk/wikis/Slurm#the-cp3slurmutils-package

.. _matplotlib: https://matplotlib.org

.. _SAMADhi: https://cp3.irmp.ucl.ac.be/samadhi/index.php

.. _virtual environment: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. _libtorch: https://pytorch.org/cppdocs/

.. _tensorflow-c: https://www.tensorflow.org/install/lang_c

.. |---| unicode:: U+2014
   :trim:
