# Bamboo

Bamboo is a python framework for HEP analyses. In this page, its explanation is mostly intended to be used in specific analyses done in Reconstruction Tools & Benchmarks Group under UPSG-CMS. For the original repo please check [https://gitlab.cern.ch/cp3-cms/bamboo](https://gitlab.cern.ch/cp3-cms/bamboo)

Documentation: [https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html)

### Setup instructions

The following set of commands would be sufficient to configure bamboo on lxplus.

```
mkdir bamboodev
cd bamboodev
```

Make a virtualenv:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```

Clone and install bamboo:

```
git clone -o upstream https://gitlab.cern.ch/aguzel/bamboo.git
python -m pip install ./bamboo
```

Clone and install plotIt:

```
git clone -o upstream https://github.com/cp3-llbb/plotIt.git
mkdir build-plotit
cd build-plotit
cmake -DCMAKE_INSTALL_PREFIX=$VIRTUAL_ENV ../plotIt
make -j4 install
cd -
```

For the .tex output, the plotit python module from [https://github.com/Oguz-Guzel/mplbplot/tree/py3compat](https://github.com/Oguz-Guzel/mplbplot/tree/py3compat) is required
```
git clone git@github.com:Oguz-Guzel/mplbplot.git
git checkout py3compat
pip install -e .
```

### Starting to use bamboo

You need to execute the following commands in the bamboo setup directory, every time logged in and want to work with bamboo:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
source bamboovenv/bin/activate
```

Testing bamboo quickly:

```
bambooRun -m bamboo/examples/nanozmumu.py:NanoZMuMu bamboo/examples/test1.yml -o test1
```

If this command works well, you can continue to start with your analysis.

### Analysis module

The analysis module of the bamboo is a python file. Please see: [/rtb/phaseII-analysis.py](https://gitlab.cern.ch/aguzel/bamboo/-/blob/master/rtb/phaseII-analysis.py).
In the analysis module, the physics objects are retreived from ROOT samples and defined via simple cuts. The selections are applied and finally the plots are defined. One can also add a cutflow report on demand which a report on selection efficiencies are printed out on terminal. 

### Analysis configuration file (YAML)

A YAML file is needed to specify the samples which are to be used in the analysis. In a configuration file, the luminosity, x-sec, sample groups and plot configurations (e.g. stacked plots, shape comparison etc.) are defined. Please check [documentation](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html) for further detail.s

### Doing analysis

The usual command line for doing analysis is as following:

```
bambooRun -m path/to/my-analsis.py path/to/configuration-file.yml -o path/to/output-directory

```

### Outputs

See examples: [http://aguzel.web.cern.ch/aguzel/rtb/](http://aguzel.web.cern.ch/aguzel/rtb/)