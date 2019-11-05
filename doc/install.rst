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

   source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
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
- the ``makePUReWeightJSON.py`` script can make a plot of the PU distributions
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

   git clone -o upstream git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git /path/to/your/bambooclone
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

Putting the above commands together, the following should give you a virtual
environment with bamboo_, and a clone of bamboo_ and plotIt in case you need to
modify them (they can be updated with ``git pull`` and ``pip install --upgrade``):

.. code-block:: sh

   mkdir bamboodev
   cd bamboodev
   # make a virtualenv
   source /cvmfs/sft.cern.ch/lcg/views/LCG_95apython3/x86_64-centos7-gcc8-opt/setup.sh
   python -m venv bamboovenv
   source bamboovenv/bin/activate
   # clone and install bamboo
   git clone -o upstream git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git
   pip install ./bamboo
   # clone and install plotIt
   git clone -o upstream https://github.com/cp3-llbb/plotIt.git
   mkdir build-plotit
   cd build-plotit
   cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
   make -j2 install
   cd -

Test your setup
---------------

TODO insert tutorial here (merge with existing)

Now you can run a few simple tests on a CMS NanoAOD to see if the installation
was successful.
First, we can pretend we are a 'worker' task, which processes trees and outputs
a file with histograms, with a test module like :py:mod:`examples.nanozmumu`:

.. code-block:: sh

   bambooRun -m /path/to/your/bambooclone/examples/nanozmumu.py:NanoZMuMu --distributed=worker --sample=DY_M50 /path/to/your/bambooclone/tests/data/DY_M50_2016.root -o testh1.root

(``--distributed=worker`` is needed to interpret the positional arguments as
input file names, in sequential mode (no ``--distributed`` option) and for
the driver task (``--distributed=driver``) the positional argument is reserved
for a yaml file that contains more information, such as input file locations
for several samples, normalisation etc. - there are a few examples).

The normal way of running bamboo is with an ``analysis.yml`` file:

.. code-block:: sh

   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu bamboo/examples/test_nanozmm.yml --envConfig=bamboo/examples/ingrid.ini -o test_nanozmm_1

If all went well, you should have a (very low-statistics) dimuon Z peak plot in
``test_nanozmm_1/plots/dimu_M.pdf``.
There is also an example that does the same with about 3.1/fb of CMS dimuon data
from 2016 (please note that this needs a valid grid proxy to retrieve the file
lists, and the files to be avialable locally under the ``storageroot`` directory
of the ``[das]`` section in the argument to ``-envConfig``):

.. code-block:: sh

   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu bamboo/examples/analysis_zmm.yml --envConfig=bamboo/examples/ingrid.ini -o test_nanozmm_2

To run this on slurm it is sufficient to add ``--distributed=driver`` (a task
with two jobs will be created, one for each sample sample).


.. _bamboo: https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html

.. _CP3SlurmUtils: https://cp3-git.irmp.ucl.ac.be/cp3-support/helpdesk/wikis/Slurm#the-cp3slurmutils-package

.. _matplotlib: https://matplotlib.org

.. _SAMADhi: https://cp3.irmp.ucl.ac.be/samadhi/index.php

.. _virtual environment: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. |---| unicode:: U+2014
   :trim:
