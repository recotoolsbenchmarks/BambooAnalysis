"""
Analysis helper functions that don't fit in treefunctions
"""
import copy
import logging
logger = logging.getLogger(__name__)
import os.path
import subprocess
import urllib.parse
import yaml

def addLumiMask(sel, jsonName, runRange=None, runAndLS=None, name="goodlumis"):
    from . import treefunctions as op
    """ Refine selection with a luminosity block filter

    Typically applied directly to the root selection (for data).
    runAndLS should be a tuple of expressions with the run number and luminosity block ID.
    The run range is used to limit the part of the JSON file to consider,
    see the LumiMask helper class for details.
    """
    if runAndLS is None:
        raise RuntimeError("Cannot construct a filter for the good lumi blocks without accessors (backend.create(..., runAndLS=XXX)), tree->(run, LS)")
    lumiSel = op.define("LumiMask", 'const auto <<name>> = LumiMask::fromJSON("{0}"{1});'.format(
                jsonName, (", {0:d}, {1:d}".format(*runRange) if runRange is not None else "")))
    return sel.refine(name, cut=lumiSel.accept(*runAndLS))

def downloadCertifiedLumiFiles(taskArgs):
    """ download certified lumi files (if needed) and replace in args """
    taskArgs = copy.deepcopy(taskArgs)
    certifLumiFiles = set(kwargs["certifiedLumiFile"] for args,kwargs in taskArgs)
    ## download if needed
    clf_downloaded = dict()
    for clfu in certifLumiFiles:
        purl = urllib.parse.urlparse(clfu)
        if purl.scheme in ("http", "https"):
            fname = purl.path.split("/")[-1]
            if os.path.exists(fname):
                logger.warning("File {0} exists, it will not be downloaded again from {1}".format(fname, clfu))
            else:
                subprocess.check_call(["wget", clfu], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            clf_downloaded[clfu] = fname
    ## update args
    for args,kwargs in taskArgs:
        if "certifiedLumiFile" in kwargs:
            clf = kwargs["certifiedLumiFile"]
            if clf in clf_downloaded:
                kwargs["certifiedLumiFile"] = clf_downloaded[clf]

    return taskArgs, set(clf_downloaded.keys())

def parseAnalysisConfig(anaCfgName, redodbqueries=False, overwritesamplefilelists=False, envConfig=None):
    cfgDir = os.path.dirname(os.path.abspath(anaCfgName))
    with open(anaCfgName) as anaCfgF:
        analysisCfg = yaml.load(anaCfgF)
    ## finish loading samples (file lists)
    samples = dict()
    for smpName, smpCfg in analysisCfg["samples"].items():
        smp = copy.deepcopy(smpCfg)
        ## read cache, if it's there
        listfile, cachelist = None, []
        if "files" in smpCfg and str(smpCfg["files"]) == smpCfg["files"]:
            listfile = smpCfg["files"] if os.path.isabs(smpCfg["files"]) else os.path.join(cfgDir, smpCfg["files"])
            if os.path.isfile(listfile):
                with open(listfile) as smpF:
                    cachelist = [ fn for fn in [ ln.strip() for ln in smpF ] if len(fn) > 0 ]

        if "db" in smpCfg and ( "files" not in smpCfg or len(cachelist) == 0 or redodbqueries ):
            if ":" not in smpCfg["db"]:
                raise RuntimeError("'db' entry should be of the format 'protocol:location', e.g. 'das:/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'")
            protocol, dbLoc = smpCfg["db"].split(":")
            files = []
            if protocol == "das":
                dasConfig = envConfig["das"]
                dasQuery = "file dataset={0}".format(dbLoc)
                files = [ os.path.join(dasConfig["storageroot"], fn) for fn in [ ln.strip() for ln in subprocess.check_output(["dasgoclient", "-query", dasQuery]).split() ] if len(fn) > 0 ]
                if len(files) == 0:
                    raise RuntimeError("No files found with DAS query {0}".format(dasQuery))
                ## TODO improve: check that files are locally available, possibly fall back to xrootd otherwise; check for grid proxy before querying; maybe do queries in parallel
            elif protocol == "samadhi":
                logger.warning("SAMADhi queries are not implemented yet")
            else:
                raise RuntimeError("Unsupported protocol in '{0}': {1}".format(smpCfg["db"], protocol))
            smp["files"] = files
            if listfile and ( len(cachelist) == 0 or overwritesamplefilelists ):
                with open(listfile, "w") as listF:
                    listF.writelines(files)
        elif "files" not in smpCfg:
            raise RuntimeError("Cannot load files for {0}: neither 'db' nor 'files' specified".format(smpName))
        elif listfile:
            if len(cachelist) == 0:
                raise RuntimeError("No file names read from {0}".format())
            smp["files"] = cachelist
        else: ## list in yml
            smp["files"] = smpCfg["files"]
        samples[smpName] = smp
    analysisCfg["samples"] = samples
    return analysisCfg

def readEnvConfig(explName=None):
    """ Read computing environment config file (batch system, storage site etc.)

    For using a batch cluster, the [batch] section should have a 'backend' key,
    and there should be a section with the name of the backend (slurm, htcondor...),
    see bamboo.batch_<backend> for details.
    The storage site information needed to resolve the PFNs for datasets retrieved from DAS
    should be specified under the [das] section (sitename and storageroot).
    """
    import os
    from configparser import ConfigParser
    def readFromFile(name):
        cfgp = ConfigParser()
        cfgp.read(name)
        cfg = dict((sName, dict(cfgp[sName])) for sName in cfgp.sections())
        return cfg

    xdgCfg = os.getenv("XDG_CONFIG_HOME", os.path.join(os.path.expanduser("~"), ".config"))
    toTry = ["bamboo.ini", "bamboorc", os.path.join(xdgCfg, "bamboorc")]
    if explName:
        toTry.insert(0, explName)
    for iniName in toTry:
        if os.path.exists(iniName):
            try:
                res = readFromFile(iniName)
                logger.info("Read config from file {0}".format(iniName))
                return res
            except Exception as ex:
                logger.warning("Problem reading config file {0}: {1}".format(iniName, ex))
    raise RuntimeError("No valid config file found")

plotit_plotdefaults = {
        "x-axis"           : lambda p : "{0}".format(p.axisTitles[0]),
        "y-axis"           : "Evt",
        "y-axis-format"    : "%1% / %2$.0f",
        "normalized"       : False,
        "x-axis-format"    : lambda p : [p.binnings[0].minimum, p.binnings[0].maximum],
        "log-y"            : "both",
        "y-axis-show-zero" : True,
        "save-extensions"  : ["pdf"],
        "show-ratio"       : True,
        "sort-by-yields"   : False,
        }
def runPlotIt(config, plotList, workdir=".", resultsdir=".", plotIt="plotIt", plotDefaults=None):
    plotitCfg = (copy.deepcopy(config["plotIt"]) if "plotIt" in config else dict())
    plotitCfg["configuration"]["root"] = resultsdir
    plotit_files = dict()
    for smpN, smpCfg in config["samples"].items():
        smpOpts = dict()
        smpOpts["group"] = smpCfg["group"]
        isMC = ( smpCfg["group"] != "data" )
        smpOpts["type"] = ("mc" if isMC else "data")
        if isMC:
            smpOpts["cross-section"] = smpCfg["cross-section"]
            ## TODO add/get "generated-events" from sum of weights (can also store in the file and read with plotIt?)
        plotit_files["{0}.root".format(smpN)] = smpOpts
    plotitCfg["files"] = plotit_files
    plotit_plots = dict()
    for plot in plotList:
        plotOpts = dict(plotit_plotdefaults)
        if plotDefaults:
            plotOpts.update(plotDefaults)
        plotOpts.update(plot.plotopts)
        plotOpts = dict((k, (v(plot) if hasattr(v, "__call__") else v)) for k,v in plotOpts.items())
        plotit_plots[plot.name] = plotOpts
    plotitCfg["plots"] = plotit_plots
    cfgName = os.path.join(workdir, "plots.yml")
    with open(cfgName, "w") as plotitFile:
        yaml.dump(plotitCfg, plotitFile)

    plotsdir = os.path.join(workdir, "plots")
    if os.path.exists(plotsdir):
        logger.warning("Directory '{0}' already exists, previous plots will be overwritten".format(plotsdir))
    else:
        os.makedirs(plotsdir)

    try:
        with open(os.path.join(plotsdir, "out.log"), "w") as logFile:
            subprocess.check_call([plotIt, "-i", workdir, "-o", plotsdir, cfgName], stdout=logFile)
        logger.info("plotIt output is available in {0}".format(plotsdir))
    except subprocess.CalledProcessError as ex:
        logger.error("Command '{0}' failed with exit code {1}\n{2}".format(" ".join(ex.cmd), ex.returncode, ex.output))
