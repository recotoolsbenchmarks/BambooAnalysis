# Bamboo

Bamboo is a python framework for HEP analyses. In this page, its explanation is mostly intended to be used in specific analyses done in Reconstruction Tools & Benchmarks Group under UPSG-CMS. For the original repo please check [https://gitlab.cern.ch/cp3-cms/bamboo](https://gitlab.cern.ch/cp3-cms/bamboo)

Documentation: [https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html)

### Setup instructions

The following set of commands would be sufficient to configure bamboo on lxplus.

```
mkdir bamboo
cd bamboo
```

Make a virtualenv:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```

Clone and install bamboo:

```
git clone -o upstream git@github.com:recotoolsbenchmarks/BambooAnalysis.git
pip install ./BambooAnalysis
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

For the .tex output, the plotit python module from [https://github.com/pieterdavid/mplbplot/tree/py3compat](https://github.com/pieterdavid/mplbplot/tree/py3compat) is required
```
git clone git@github.com:pieterdavid/mplbplot.git
cd mplbplot
git checkout py3compat
pip install -e .
cd -
```

### Starting to use bamboo

You need to execute the following commands in the bamboo setup directory, every time logged in and want to work with bamboo:

```
source /cvmfs/sft.cern.ch/lcg/views/LCG_97python3/x86_64-centos7-gcc9-opt/setup.sh
python -m venv bamboovenv
source bamboovenv/bin/activate
```

### Analysis module

The analysis module of the bamboo is a python file. Please see: [/rtb/phaseII-analysis.py](https://github.com/Oguz-Guzel/bamboo/-/blob/master/rtb/phaseII-analysis.py).
In the analysis module, the physics objects are retreived from ROOT trees and defined via simple cuts. The selections are applied and finally the plots are defined. One can also add a cutflow report on demand which prints selection efficiencies in a .tex file as well as in the terminal screen. This cutflow report comes along with a histogram in .png format available in the working directory.

### Analysis configuration file (YAML)

A YAML file is needed to specify the samples to be used in the analysis. In a configuration file, the luminosity, x-sec, sample groups and plot configurations (e.g. stacked plots, shape comparison etc.) are defined. Please check [documentation](https://cp3.irmp.ucl.ac.be/~pdavid/bamboo/index.html) for further details. See an example in the /rtb.

### RTB analysis

DelphesFlat
```
bambooRun -m BambooAnalysis/rtb/phaseII-analysis.py:CMSPhase2SimTest BambooAnalysis/rtb/phaseII-analysis-FS-Delphes.yml -o output-Delphes
```
Fullsim
```
bambooRun -m BambooAnalysis/rtb/phaseII-analysis.py:CMSPhase2SimTest BambooAnalysis/rtb/phaseII-analysis-FS.yml -o output-FS
```
Delphes vs. FS
```
bambooRun -m BambooAnalysis/rtb/phaseII-analysis.py:CMSPhase2SimTest BambooAnalysis/rtb/phaseII-analysis-FS-Delphes.yml -o output-FS-Delphes
```

The pdf version of the cutflow report is obtained via:
```
cd output-FS-Delphes
pdflatex yields_HL-LHC.tex
```

### Outputs

See examples: [http://aguzel.web.cern.ch/aguzel/rtb/](http://aguzel.web.cern.ch/aguzel/rtb/)
