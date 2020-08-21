"""
The :py:mod:`bamboo.root` module collects a set of thin wrappers around ROOT
methods, and centralizes the import of the Cling interpreter global namespace
in PyROOT. For compatibility, it is recommended that user code uses
``from bamboo.root import gbl`` rather than ``import ROOT as gbl`` or
``from cppyy import gbl``.
"""
import pkg_resources

## workaround for https://sft.its.cern.ch/jira/browse/ROOT-10304, which affects only ROOT 6.18/XX
import subprocess
_rootVersion = subprocess.check_output(["root-config", "--version"]).decode().strip()
_rootVersion_split = tuple((int(iv) for tk in _rootVersion.split("/") for iv in tk.split(".")))
if _rootVersion_split[0] == 6 and _rootVersion_split[1] == 18:
    import ROOT as gbl
    gbl.PyConfig.IgnoreCommandLineOptions = True
else:
    from cppyy import gbl

def addIncludePath(incPath):
    """ Add an include path to the ROOT interpreter """
    gbl.gInterpreter.AddIncludePath(incPath)
def loadHeader(headerName):
    """ Include a C++ header in the ROOT interpreter """
    gbl.gROOT.ProcessLine('#include "{}"'.format(headerName))
def addDynamicPath(libPath):
    """ Add a dynamic library path to the ROOT interpreter"""
    gbl.gSystem.AddDynamicPath(libPath)
def loadLibrary(libName):
    """ Load a shared library in the ROOT interpreter """
    st = gbl.gSystem.Load(libName)
    if st == -1:
        raise RuntimeError("Library {0} could not be found".format(libName))
    elif st == -2:
        raise RuntimeError("Version match for library {0}".format(libName))
def findLibrary(libName):
    """ Check if a library can be found, and returns the path in that case """
    ret = gbl.gSystem.FindDynamicLibrary(gbl.TString(libName), True)
    if ret != "":
        return str(ret)

import functools
class once:
    """ Function decorator to make sure things are not loaded more than once """
    def __init__(self, fun):
        self.fun = fun
        self.has_run = False
        functools.update_wrapper(self, fun)
    def __call__(self):
        if not self.has_run:
            self.has_run = True
            return self.fun()

@once
def loadBambooExtensions():
    # Add extension libraries and necessary header files to the ROOT interpreter
    import sys
    import os.path
    pkgRoot = pkg_resources.get_distribution("bamboo").location
    instInclude = os.path.join(pkgRoot, "bamboo", "include")
    if os.path.isdir(instInclude): ## installed mode
        addIncludePath(instInclude)
        libDir = pkgRoot
    else: ## non-installed mode
        libDir = os.path.join(pkgRoot, "build", "lib")
        if not os.path.isdir(libDir):
            raise RuntimeError("No directory {0} so running in local mode, but then build/lib need to be present. Did you run 'python setup.py build'?".format(libDir))
        addIncludePath(os.path.join(pkgRoot, "build", "include"))
        addIncludePath(os.path.join(pkgRoot, "cpp"))
    addDynamicPath(libDir)
    ## now load default headers and libraries
    loadHeader("Math/VectorUtil.h")
    loadLibrary("libBambooLumiMask")
    loadLibrary("libBambooRandom")
    loadLibrary("libBinnedValues")
    for fname in ("bamboohelpers.h", "range.h", "LumiMask.h", "bamboorandom.h", "scalefactors.h", "BTagCalibrationStandalone.h"):
        loadHeader(fname)

@once
def loadJMESystematicsCalculators():
    loadBambooExtensions()
    loadLibrary("libJMEObjects")
    loadHeader("JMESystematicsCalculators.h")
    getattr(gbl, "JetVariationsCalculator::result_t") ## trigger dictionary generation

@once
def loadRochesterCorrectionCalculator():
    loadBambooExtensions()
    loadLibrary("libRoccoR")
    loadHeader("RochesterCorrectionCalculator.h")
    getattr(gbl, "RochesterCorrectionCalculator::result_t") ## trigger dictionary generation

@once
def loadLibTorch():
    if loadTensorflowC.has_run:
        raise RuntimeError("If loading both Tensorflow-C and libtorch, the latter should be loaded first")
    loadBambooExtensions()
    torchLib = findLibrary("libtorch")
    if ( not torchLib ) and pkg_resources.resource_filename("torch", "lib"):
        addDynamicPath(pkg_resources.resource_filename("torch", "lib"))
    loadLibrary("libBambooTorch")
    loadHeader("bambootorch.h")

@once
def loadTensorflowC():
    loadBambooExtensions()
    tfLib = findLibrary("libtensorflow")
    import os
    if ( not tfLib ) and "VIRTUAL_ENV" in os.environ:
        addDynamicPath(os.path.join(os.environ["VIRTUAL_ENV"], "lib"))
    loadLibrary("libtensorflow")
    loadLibrary("libBambooTensorflowC")
    loadHeader("bambootensorflowc.h")

@once
def loadlwtnn():
    loadBambooExtensions()
    lwtnnLib = findLibrary("liblwtnn")
    import os
    if ( not lwtnnLib ) and "VIRTUAL_ENV" in os.environ:
        addDynamicPath(os.path.join(os.environ["VIRTUAL_ENV"], "lib"))
    loadLibrary("liblwtnn")
    loadLibrary("libBambooLwtnn")
    loadHeader("bamboolwtnn.h")
