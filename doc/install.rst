Installation and setup
======================

Bamboo_ only depends on python3 (with pip/setuptools to install PyYAML if needed)
and a recent version of ROOT (at least 6.14, for the necessary features of
RDataFrame).
On ingrid and lxplus (or any machine with cvmfs), an easy way to get such
a recent version of ROOT is through a CMSSW release that depends on it (``10_4``
has a recent enough ROOT, but no python3 yet), or from the lcgsoft distribution,
e.g.

.. code-block:: sh

   source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh
   python -m venv bamboovenv
   source bamboovenv/bin/activate

(the second command creates a `virtual environment`_
to install python packages in, after installation it is sufficient to run two
other commands, to pick up the correct base system and then the installed
packages).

.. note::

   Querying DAS relies on dasgoclient (and a valid grid proxy). If you use this
   feature the safest is currently to run ``cms_env`` before the commands above
   (together with ``voms-proxy-init`` |---| which can alternatively also be run
   from a different shell on the same machine).

Bamboo_ can be installed in the `virtual environment`_ with pip, minimally either
one of

.. code-block:: sh

   pip install git+https://gitlab.cern.ch/cp3-cms/bamboo.git
   pip install git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git

but in the current development stage it may be useful to install from
a local clone, such that you can use it to test and propose changes, using

.. code-block:: sh

   git clone -o upstream git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git /path/to/your/bambooclone
   # copy and patch some jet-related classes from CMSSW (requires cvmfs) --- temporary
   pushd /path/to/your/bambooclone/ext
   ./getjetclasses.sh
   popd
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
which is used in many analyses, the plotIt_
tool is used. It can be installed with

.. code-block:: sh

   git clone -o upstream https://github.com/cp3-llbb/plotIt.git /path/to/your/plotitclone
   cd /path/to/your/plotitclone
   cd external
   ./build-external.sh
   cd -
   BOOST_ROOT=$CMAKE_PREFIX_PATH make -j4
   cp plotIt bamboovenv/bin

where the last command copies the ``plotIt`` executable inside the virtualenv
executable directory such that it is picked up automatically (alternatively, its
path can be passed to ``bambooRun`` with the ``--plotIt`` command-line option).

Putting the above commands together, the following should give you a virtual
environment with bamboo_, and a clone of bamboo_ and plotIt in case you need to
modify them (they can be updated with ``git pull`` and ``pip install --upgrade``):

.. code-block:: sh

   mkdir bamboodev
   cd bamboodev
   source /cvmfs/sft.cern.ch/lcg/views/LCG_94python3/x86_64-centos7-gcc8-opt/setup.sh
   python -m venv bamboovenv
   source bamboovenv/bin/activate
   git clone -o upstream git+ssh://git@gitlab.cern.ch:7999/cp3-cms/bamboo.git
   pip install ./bamboo
   git clone -o upstream https://github.com/cp3-llbb/plotIt.git
   cd plotIt/external
   ./build-external.sh
   cd ..
   BOOST_ROOT=$CMAKE_PREFIX_PATH make -j4
   cd ..
   cp plotIt/plotIt bamboovenv/bin

Now you can run a few simple tests on a CMS NanoAOD (on ingrid you could use
``/home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root``)
to see if the installation was successful.
First, we can pretend we are a 'worker' task, which processes trees and outputs
a file with histograms, with a test module like :py:mod:`examples.nanozmumu`:

.. code-block:: sh

   bambooRun -m /path/to/your/clone/examples/nanozmumu.py:NanoZMuMu --distributed=worker /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root -o testh1.root

(``--distributed=worker`` is needed to interpret the positional arguments as
input file names, in sequential mode (no ``--distributed`` option) and for
the driver task (``--distributed=driver``) the positional argument is reserved
for a json/yaml file that contains more information, such as input file
locations for several samples, normalisation etc. - there are a few examples).

A more complete example would run from an ``analysis.yml`` file (copy it to
``bamboo/examples`` because ``test_nanozmm1.yml`` specifies it as a local file
with relative path):

.. code-block:: sh

   cp /home/ucl/cp3/pdavid/bamboodev/bamboo/examples/NanoAOD_SingleMu_test.root bamboo/examples
   bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu bamboo/examples/test_nanozmm1.yml --envConfig=bamboo/examples/ingrid.ini -o test_nanozmm1

if all went well, you should have a dimuon Z peak plot in
``test_nanozmm1/plots/dimu_M.pdf``. To run on slurm add
``--distributed=driver``.

Passing the ``--envConfig`` option can in practice be avoided by copying the
appropriate file to ``~/.config/bamboorc``. It is necessary to pick up the
configuration of the computing environment (files for ingrid and lxplus are
included), e.g. how to access the file storage, which batch submission system
to use (currently slurm and HTCondor are supported). Bamboo_ tries to find it
first from the ``--envConfig`` option, then from ``bamboo.ini`` in the current
directory, then ``$XDG_CONFIG_HOME/bamboorc`` (which typically resolves to
``~/.config/bamboorc``).

.. _bamboo: http://to-fill-bamboodocs-home

.. _virtual environment: https://packaging.python.org/tutorials/installing-packages/#creating-virtual-environments

.. _plotIt: https://github.com/cp3-llbb/plotIt

.. |---| unicode:: U+2014
   :trim:
