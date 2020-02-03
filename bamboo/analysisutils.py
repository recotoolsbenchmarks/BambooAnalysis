""" The :py:mod:`bamboo.analysisutils` module bundles a number of more
specific helper methods that use the tree decorators and integrate with
other components, connect to external services, or are factored out of the
classes in :py:mod:`bamboo.analysismodules` to facilitate reuse.
"""
import copy
import logging
logger = logging.getLogger(__name__)
import os.path
import subprocess
import urllib.parse
import yaml

_SAMADhi_found = False
try:
    from cp3_llbb.SAMADhi.SAMADhi import Sample, SAMADhiDB
    _SAMADhi_found = True
except ImportError as ex:
    logger.warning("Could not load SAMADhi, please install the SAMADhi library if you want to use the database to locate samples")

bamboo_cachedir = os.path.join(os.getenv("XDG_CACHE_HOME", os.path.join(os.path.expanduser("~"), ".cache")), "bamboo")

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

def downloadCertifiedLumiFiles(taskArgs, workdir="."):
    """ download certified lumi files (if needed) and replace in args """
    taskArgs = copy.deepcopy(taskArgs)
    certifLumiFiles = set(kwargs["certifiedLumiFile"] for args,kwargs in taskArgs if "certifiedLumiFile" in kwargs)
    ## download if needed
    clf_downloaded = dict()
    for clfu in certifLumiFiles:
        purl = urllib.parse.urlparse(clfu)
        if purl.scheme in ("http", "https"):
            fname = os.path.join(workdir, purl.path.split("/")[-1])
            if os.path.exists(fname):
                logger.warning("File {0} exists, it will not be downloaded again from {1}".format(fname, clfu))
            else:
                subprocess.check_call(["wget", "--directory-prefix={0}".format(workdir), clfu], stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
            clf_downloaded[clfu] = os.path.abspath(fname)
    ## update args
    for args,kwargs in taskArgs:
        if "certifiedLumiFile" in kwargs:
            clf = kwargs["certifiedLumiFile"]
            if clf in clf_downloaded:
                kwargs["certifiedLumiFile"] = clf_downloaded[clf]

    return taskArgs, set(clf_downloaded.keys())

def _dasLFNtoPFN(lfn, dasConfig):
    localPFN = os.path.join(dasConfig["storageroot"], lfn.lstrip("/"))
    if dasConfig.get("checklocalfiles", False) and "xrootdredirector" in dasConfig:
        if os.path.isfile(localPFN):
            return localPFN
        else:
            xrootdPFN = "root://{redirector}//{lfn}".format(redirector=dasConfig["xrootdredirector"], lfn=lfn)
            logger.warning("PFN {0} not available, falling back to xrootd with {1}".format(localPFN, xrootdPFN))
            return xrootdPFN
    else:
        return localPFN

def sample_resolveFiles(smpCfg, redodbqueries=False, overwritesamplefilelists=False, envConfig=None, cfgDir="."):
    smp = copy.deepcopy(smpCfg)
    ## read cache, if it's there
    listfile, cachelist = None, []
    if "files" in smpCfg and str(smpCfg["files"]) == smpCfg["files"]:
        listfile = smpCfg["files"] if os.path.isabs(smpCfg["files"]) else os.path.join(cfgDir, smpCfg["files"])
        if os.path.isfile(listfile):
            with open(listfile) as smpF:
                cachelist = [ fn for fn in [ ln.strip() for ln in smpF ] if len(fn) > 0 ]

    if "db" in smpCfg and ( "files" not in smpCfg or len(cachelist) == 0 or redodbqueries ):
        files = []
        for dbEntry in (smpCfg["db"] if str(smpCfg["db"]) != smpCfg["db"] else [smpCfg["db"]]): ## convert to list if string
            if ":" not in dbEntry:
                raise RuntimeError("'db' entry should be of the format 'protocol:location', e.g. 'das:/SingleMuon/Run2016E-03Feb2017-v1/MINIAOD'")
            protocol, dbLoc = dbEntry.split(":")
            if protocol == "das":
                dasConfig = envConfig["das"]
                dasQuery = "file dataset={0}".format(dbLoc)
                entryFiles = [ _dasLFNtoPFN(lfn, dasConfig) for lfn in [ ln.strip() for ln in subprocess.check_output(["dasgoclient", "-query", dasQuery]).decode().split() ] if len(lfn) > 0 ]
                files += entryFiles
                if len(entryFiles) == 0:
                    raise RuntimeError("No files found with DAS query {0}".format(dasQuery))
                ## TODO improve: check for grid proxy before querying; maybe do queries in parallel
            elif protocol == "samadhi":
                if not _SAMADhi_found:
                    raise RuntimeError("SAMADhi could not be found, cannot resolve '{0}'".format(dbEntry))
                samaCred = "~/.samadhi"
                if "SAMADhi" in envConfig and "credentials" in envConfig["SAMADhi"]:
                    samaCred = envConfig["SAMADhi"]["credentials"]
                with SAMADhiDB(credentials=samaCred) as db:
                    if dbLoc.isnumeric():
                        descr = "id {0}".format(dbLoc)
                        sample = Sample.get_or_none(Sample.id == int(dbLoc))
                    else:
                        descr = "name '{0}'".format(dbLoc)
                        sample = Sample.get_or_none(Sample.name == dbLoc)
                    if not sample:
                        raise RuntimeError("Could not find sample with {0} in SAMADhi".format(descr))
                    entryFiles = [ f.pfn for f in sample.files ]
                    if len(entryFiles) == 0:
                        raise RuntimeError("No files found with SAMADhi {0}".format(descr))
                files += entryFiles
            else:
                raise RuntimeError("Unsupported protocol in '{0}': {1}".format(dbEntry, protocol))
        smp["files"] = files
        if listfile and ( len(cachelist) == 0 or overwritesamplefilelists ):
            with open(listfile, "w") as listF:
                listF.write("\n".join(files))
    elif "files" not in smpCfg:
        raise RuntimeError("Cannot load files for {0}: neither 'db' nor 'files' specified".format(smpName))
    elif listfile:
        if len(cachelist) == 0:
            raise RuntimeError("No file names read from {0}".format())
        smp["files"] = cachelist
    else: ## list in yml
        smp["files"] = [ (fn if os.path.isabs(fn) or urllib.parse.urlparse(fn).scheme != "" in fn else os.path.join(cfgDir, fn)) for fn in smpCfg["files"] ]
    return smp

class YMLIncludeLoader(yaml.SafeLoader):
    """Custom yaml loading to support including config files. Use `!include (file)` to insert content of `file` at that position."""
    
    def __init__(self, stream):
        super(YMLIncludeLoader, self).__init__(stream)
        self._root = os.path.split(stream.name)[0]

    def include(self, node):
        filename = os.path.join(self._root, self.construct_scalar(node))
        with open(filename, 'r') as f:
            return yaml.load(f, YMLIncludeLoader)

YMLIncludeLoader.add_constructor('!include', YMLIncludeLoader.include)

def parseAnalysisConfig(anaCfgName, resolveFiles=True, redodbqueries=False, overwritesamplefilelists=False, envConfig=None):
    cfgDir = os.path.dirname(os.path.abspath(anaCfgName))
    with open(anaCfgName) as anaCfgF:
        analysisCfg = yaml.load(anaCfgF, YMLIncludeLoader)
    if resolveFiles:
        analysisCfg["samples"] = dict((smpName,
            sample_resolveFiles(smpCfg, redodbqueries=redodbqueries, overwritesamplefilelists=overwritesamplefilelists, envConfig=envConfig, cfgDir=cfgDir))
            for smpName, smpCfg in analysisCfg["samples"].items())
    return analysisCfg

def getAFileFromAnySample(samples, redodbqueries=False, overwritesamplefilelists=False, envConfig=None, cfgDir="."):
    """ Helper method: get a file from any sample (minimizing the risk of errors)

    Tries to find any samples with:
    - a list of files
    - a cache file
    - a SAMADhi path
    - a DAS path

    If successful, a single read / query is sufficient to retrieve a file
    """
    ## list of files -> return 1st
    for smpNm,smpCfg in samples.items():
        if ( "files" in smpCfg ) and ( str(smpCfg["files"]) != smpCfg["files"] ):
            return smpNm,smpCfg,(smpCfg["files"][0] if (os.path.isabs(smpCfg["files"][0]) or urllib.parse.urlparse(smpCfg["files"][0]).scheme) else os.path.join(cfgDir, smpCfg["files"][0]))
    ## try to get them from a cache file or database (ordered by less-to-more risky)
    failed_names = set()
    for method, condition in [
            (" from cache file", (lambda smpCfg : "files" in smpCfg and ( str(smpCfg["files"]) == smpCfg["files"] ))),
            (" from SAMADhi"   , (lambda smpCfg : "db" in smpCfg and smpCfg["db"].startswith("samadhi:"))),
            (" from DAS"       , (lambda smpCfg : "db" in smpCfg and smpCfg["db"].startswith("das:"))),
            (""                , (lambda smpCfg : True))
            ]:
        for smpNm,smpCfg in samples.items():
            if smpNm not in failed_names and condition(smpCfg):
                try:
                    smpCfg = sample_resolveFiles(smpCfg, redodbqueries=redodbqueries, overwritesamplefilelists=overwritesamplefilelists, envConfig=envConfig, cfgDir=cfgDir)
                    return smpNm,smpCfg,smpCfg["files"][0]
                except Exception as ex:
                    failed_names.add(smpNm)
                    logger.warning("Problem while resolving files for {0}{1}: {2!s}".format(smpNm, method, ex))

    raise RuntimeError("Failed to resolve a file from any sample (see the warnings above for more information)")

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
    logger.error("No valid environment config file found, please copy one from the .ini files from the examples directory to ~/.config/bamboorc, or pass the --envConfig option")

plotit_plotdefaults = {
        "x-axis"           : lambda p : "{0}".format(p.axisTitles[0]),
        "x-axis-range"     : lambda p : [p.binnings[0].minimum, p.binnings[0].maximum],
        }
def runPlotIt(config, plotList, workdir=".", resultsdir=".", plotIt="plotIt", plotDefaults=None, readCounters=lambda f : -1., eras=("all", None), verbose=False):
    eraMode, eras = eras
    if not eras: ## from config if not specified
        eras = list(config["eras"].keys())
    plotitCfg = copy.deepcopy(config["plotIt"]) if plotIt in config else {"configuration":{}}
    plotitCfg["configuration"].update({
        "root" : os.path.relpath(resultsdir, workdir),
        "eras" : eras,
        "luminosity" : dict((era, config["eras"][era]["luminosity"]) for era in eras)
        })
    plotit_files = dict()
    for smpN, smpCfg in config["samples"].items():
        if smpCfg.get("era") in eras:
            resultsName = "{0}.root".format(smpN)
            smpOpts = dict(smpCfg)
            isMC = ( smpCfg.get("group") != "data" )
            if "type" not in smpOpts:
                smpOpts["type"] = ("mc" if isMC else "data")
            if isMC:
                if "cross-section" not in smpCfg:
                    logger.warning("Sample {0} is of type MC, but no cross-section specified".format(smpN))
                smpOpts["cross-section"] = smpCfg.get("cross-section", 1.)
                from .root import gbl
                resultsFile = gbl.TFile.Open(os.path.join(resultsdir, resultsName))
                if "generated-events" not in smpCfg:
                    logger.error(f"No key 'generated-events' found for MC sample {smpNm}, normalization will be wrong")
                else:
                    try:
                        smpOpts["generated-events"] = float(smpCfg["generated-events"])
                    except ValueError: ## not float
                        try:
                            counters = readCounters(resultsFile)
                            smpOpts["generated-events"] = counters[smpCfg["generated-events"]]
                        except Exception as ex:
                            logger.error("Problem reading counters for sample {0} (file {1}), normalization may be wrong (exception: {2!r})".format(smpN, os.path.join(resultsdir, resultsName), ex))
            plotit_files[resultsName] = smpOpts
    plotitCfg["files"] = plotit_files
    plotit_plots = dict()
    for plot in plotList:
        if len(plot.variables) == 1:
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

    out_extraOpts = []
    if len(eras) > 1 and eraMode in ("all", "combined"):
        out_extraOpts.append((os.path.join(workdir, "plots"), []))
    if len(eras) == 1 or eraMode in ("split", "all"):
        for era in eras:
            out_extraOpts.append((os.path.join(workdir, "plots_{0}".format(era)), ["-e", era]))
    for plotsdir, extraOpts in out_extraOpts:
        if os.path.exists(plotsdir):
            logger.warning("Directory '{0}' already exists, previous plots will be overwritten".format(plotsdir))
        else:
            os.makedirs(plotsdir)
        try:
            plotItLog = os.path.join(plotsdir, "out.log")
            plotItArgs = [plotIt, "-i", workdir, "-o", plotsdir]+extraOpts+[cfgName]
            if verbose:
                logger.debug("Running command `{0}`, with logfile {1}".format(" ".join(plotItArgs), plotItLog))
            with open(plotItLog, "w") as logFile:
                subprocess.check_call(plotItArgs, stdout=logFile)
            logger.info("plotIt output is available in {0}".format(plotsdir))
        except subprocess.CalledProcessError as ex:
            logger.error("Command '{0}' failed with exit code {1}\n{2}".format(" ".join(ex.cmd), ex.returncode, ex.output))

def configureJets(variProxy, jetType, jec=None, jecLevels="default", smear=None, useGenMatch=True, genMatchDR=0.2, genMatchDPt=3., jesUncertaintySources=None, cachedir=None, mayWriteCache=False, enableSystematics=None, isMC=False, backend=None, uName=""):
    """ Reapply JEC, set up jet smearing, or prepare JER/JES uncertainties collections

    :param variProxy: jet variations proxy, e.g. ``tree._Jet``
    :param jetType: jet type, e.g. AK4PFchs
    :param smear: tag of resolution (and scalefactors) to use for smearing (no smearing is done if unspecified)
    :param jec: tag of the new JEC to apply, or for the JES uncertainties (pass an empty list to jecLevels to produce only the latter without reapplying the JEC)
    :param jesUncertaintySources: list of jet energy scale uncertainty sources (see `the JECUncertaintySources twiki page <https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources>`_)
    :param enableSystematics: systematics variations to enable (default: all that are available; pass an empty set to disable)

    :param useGenMatch: use matching to generator-level jets for resolution smearing
    :param genMatchDR: DeltaR for generator-level jet matching (half the cone size is recommended, default is 0.2)
    :param genMatchDPt: maximal relative PT difference (in units of the resolution) between reco and gen jet
    :param jecLevels: list of JEC levels to apply (if left out the recommendations are used: L1FastJet, L2Relative, L3Absolute, and also L2L3Residual for data)
    :param cachedir: alternative root directory to use for the txt files cache, instead of ``$XDG_CACHE_HOME/bamboo`` (usually ``~/.cache/bamboo``)
    :param mayWriteCache: flag to indicate if this task is allowed to write to the cache status file (set to False for worker tasks to avoid corruption due to concurrent writes)

    :param isMC: MC or not
    :param backend: backend pointer (returned from :py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree`)
    :param uName: unique name for the correction calculator (sample name is a safe choice)
    """
    from . import treefunctions as op # first to load default headers/libs, if still needed
    from .treefunctions import _to, _tp ## treeoperations and treeproxies
    from itertools import repeat
    aJet = variProxy.orig[0]
    args = [ getattr(aJet, comp).op.arg for comp in ("pt", "eta", "phi", "mass", "rawFactor", "area") ]
    args.append(_to.GetColumn("Float_t", "fixedGridRhoFastjetAll"))
    if isMC:
        evt = variProxy._parent
        args.append((evt.run<<20) + (evt.luminosityBlock<<10) + evt.event + 1 + op.static_cast("unsigned",
                    op.switch(op.rng_len(variProxy.orig) != 0, variProxy.orig[0].eta/.01, op.c_float(0.))))
        aGJet = evt.GenJet[0]
        args += [ getattr(aGJet, comp).op.arg for comp in ("pt", "eta", "phi", "mass") ]
    else:
        args.append(_tp.makeConst(0, "unsigned")) # no seed
        args += list(repeat(_to.ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}"), 4))
    ## load necessary library and header(s)
    from .root import loadJMESystematicsCalculators, gbl
    loadJMESystematicsCalculators()
    ## define calculator and initialize
    jetcalcName = backend.symbol("JetVariationsCalculator <<name>>{{}}; // for {0}".format(uName), nameHint="bamboo_jetVarCalc{0}".format("".join(c for c in uName if c.isalnum())))
    variProxy._initCalc(op.extVar("JetVariationsCalculator", jetcalcName), calcHandle=getattr(gbl, jetcalcName), args=args)
    calc = variProxy.calc
    from .jetdatabasecache import JetDatabaseCache
    if smear is not None:
        jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase", cachedir=cachedir, mayWrite=mayWriteCache)
        mcPTRes = jrDBCache.getPayload(smear, "PtResolution", jetType)
        mcResSF = jrDBCache.getPayload(smear, "SF", jetType)
        calc.setSmearing(mcPTRes, mcResSF, useGenMatch, genMatchDR, genMatchDPt)
    if jec is not None:
        if jecLevels == "default":
            # "L3Absolute" left out because it is dummy according to https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
            if jec.endswith("_DATA"):
                jecLevels = ["L1FastJet", "L2Relative", "L2L3Residual"]
            elif jec.endswith("_MC"):
                # "L2L3Residual" could be added, but it is dummy for MC according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorApplication
                jecLevels = ["L1FastJet", "L2Relative"]
            else:
                raise ValueError("JEC tag {0} does not end with '_DATA' or '_MC', so the levels cannot be guessed. Please specify the JEC levels explicitly")
        jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase", cachedir=cachedir, mayWrite=mayWriteCache)
        from .root import gbl
        if jecLevels:
            jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
            for jLev in jecLevels:
                plf = jecDBCache.getPayload(jec, jLev, jetType)
                params = gbl.JetCorrectorParameters(plf)
                jecParams.push_back(params)
            calc.setJEC(jecParams)
        if jesUncertaintySources:
            plf = jecDBCache.getPayload(jec, "UncertaintySources", jetType)
            for src in jesUncertaintySources:
                params = gbl.JetCorrectorParameters(plf, src)
                calc.addJESUncertainty(src, params)
    variProxy._initFromCalc()
    if enableSystematics is not None:
        if str(enableSystematics) == enableSystematics:
            enableSystematics = [enableSystematics]
        avail = list(variProxy.calc.available())
        enable = [ vari for vari in avail if vari != "nominal" and any(vari.startswith(ena) for ena in enableSystematics) ]
        if enable:
            for opWithSyst in variProxy.brMapMap[variProxy.withSystName].values():
                opWithSyst.variations = enable ## fine, just (re)creataed by _initFromCalc
            logger.info("Enabled systematic variations for {0} collection: {1}".format(jetsName, " ".join(enable)))

def configureType1MET(variProxy, jec=None, smear=None, useGenMatch=True, genMatchDR=0.2, genMatchDPt=3., jesUncertaintySources=None, cachedir=None, mayWriteCache=False, enableSystematics=None, isMC=False, backend=None, uName=""):
    """ Reapply JEC, set up jet smearing, or prepare JER/JES uncertainties collections

    :param variProxy: MET variations proxy, e.g. ``tree._MET``
    :param smear: tag of resolution (and scalefactors) to use for smearing (no smearing is done if unspecified)
    :param jec: tag of the new JEC to apply, or for the JES uncertainties
    :param jesUncertaintySources: list of jet energy scale uncertainty sources (see `the JECUncertaintySources twiki page <https://twiki.cern.ch/twiki/bin/viewauth/CMS/JECUncertaintySources>`_)
    :param enableSystematics: systematics variations to enable (default: all that are available; pass an empty set to disable)

    :param useGenMatch: use matching to generator-level jets for resolution smearing
    :param genMatchDR: DeltaR for generator-level jet matching (half the cone size is recommended, default is 0.2)
    :param genMatchDPt: maximal relative PT difference (in units of the resolution) between reco and gen jet
    :param cachedir: alternative root directory to use for the txt files cache, instead of ``$XDG_CACHE_HOME/bamboo`` (usually ``~/.cache/bamboo``)
    :param mayWriteCache: flag to indicate if this task is allowed to write to the cache status file (set to False for worker tasks to avoid corruption due to concurrent writes)

    :param be: backend pointer
    :param uName: unique name (sample name is a safe choice)
    :param isMC: MC or not
    """
    isFixEE2017 = variProxy.orig.pt.op.name.startswith("METFixEE2017_")
    from . import treefunctions as op # first to load default headers/libs, if still needed
    from .treefunctions import _to, _tp ## treeoperations and treeproxies
    from itertools import repeat
    evt = variProxy._parent
    jets = evt._Jet.orig if hasattr(evt, "_Jet") else evt.Jet[0]
    args = [ getattr(jets[0], comp).op.arg for comp in ("pt", "eta", "phi", "mass", "rawFactor", "area", "muonSubtrFactor", "neEmEF", "chEmEF") ]
    args.append(_to.GetColumn("Float_t", "fixedGridRhoFastjetAll"))
    evt = variProxy._parent
    if isMC:
        args.append((evt.run<<20) + (evt.luminosityBlock<<10) + evt.event + 1 + op.static_cast("unsigned",
                    op.switch(op.rng_len(jets) != 0, jets[0].eta/.01, op.c_float(0.))))
        aGJet = evt.GenJet[0]
        args += [ getattr(aGJet, comp).op.arg for comp in ("pt", "eta", "phi", "mass") ]
    else:
        args.append(_tp.makeConst(0, "unsigned")) # no seed
        args += list(repeat(_to.ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}"), 4))
    args += [ evt.RawMET.pt, evt.RawMET.pt, variProxy.orig.MetUnclustEnUpDeltaX, variProxy.orig.MetUnclustEnUpDeltaY ]
    aT1Jet = evt.CorrT1METJet[0]
    args += [ getattr(aT1Jet, comp).op.arg for comp in ("rawPt", "eta", "phi", "area", "muonSubtrFactor") ]
    args += [ getattr(aT1Jet, comp).op.arg if hasattr(aT1Jet, comp) else _to.ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}") for comp in ("neEmEF", "chEmEF") ]
    if isFixEE2017:
        args += [ evt.MET.phi, evt.MET.pt, variProxy.orig.phi, variProxy.orig.pt ]
    ## load necessary library and header(s)
    from .root import loadJMESystematicsCalculators, gbl
    loadJMESystematicsCalculators()
    ## define calculator and initialize
    calcType = "Type1METVariationsCalculator" if not isFixEE2017 else "FixEE2017Type1METVariationsCalculator"
    metcalcName = backend.symbol("{calcType} <<name>>{{}}; // for {0}".format(uName, calcType=calcType), nameHint="bamboo_Type1METVarCalc{0}".format("".join(c for c in uName if c.isalnum())))
    variProxy._initCalc(op.extVar(calcType, metcalcName), calcHandle=getattr(gbl, metcalcName), args=args)
    calc = variProxy.calc
    ## configure the calculator
    jetType = "AK4PFchs"
    from .jetdatabasecache import JetDatabaseCache
    if smear is not None:
        jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase", cachedir=cachedir, mayWrite=mayWriteCache)
        mcPTRes = jrDBCache.getPayload(smear, "PtResolution", jetType)
        mcResSF = jrDBCache.getPayload(smear, "SF", jetType)
        calc.setSmearing(mcPTRes, mcResSF, useGenMatch, genMatchDR, genMatchDPt)
    if jec is not None:
        jecLevels_L1 = ["L1FastJet"]
        # "L3Absolute" left out because it is dummy according to https://twiki.cern.ch/twiki/bin/view/CMS/IntroToJEC#Mandatory_Jet_Energy_Corrections
        if jec.endswith("_DATA"):
            jecLevels = jecLevels_L1 + ["L2Relative", "L2L3Residual"]
        elif jec.endswith("_MC"):
            # "L2L3Residual" could be added, but it is dummy for MC according to https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookJetEnergyCorrections#JetCorApplication
            jecLevels = jecLevels_L1 + ["L2Relative"]
        else:
            raise ValueError("JEC tag {0} does not end with '_DATA' or '_MC', so the levels cannot be guessed")
        jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase", cachedir=cachedir, mayWrite=mayWriteCache)
        from .root import gbl
        if jecLevels:
            jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
            for jLev in jecLevels:
                plf = jecDBCache.getPayload(jec, jLev, jetType)
                params = gbl.JetCorrectorParameters(plf)
                jecParams.push_back(params)
            calc.setJEC(jecParams)
        if jecLevels_L1:
            jecParams = getattr(gbl, "std::vector<JetCorrectorParameters>")()
            for jLev in jecLevels_L1:
                plf = jecDBCache.getPayload(jec, jLev, jetType)
                params = gbl.JetCorrectorParameters(plf)
                jecParams.push_back(params)
            calc.setL1JEC(jecParams)
        if jesUncertaintySources:
            plf = jecDBCache.getPayload(jec, "UncertaintySources", jetType)
            for src in jesUncertaintySources:
                params = gbl.JetCorrectorParameters(plf, src)
                calc.addJESUncertainty(src, params)
    variProxy._initFromCalc()
    if enableSystematics is not None:
        if str(enableSystematics) == enableSystematics:
            enableSystematics = [enableSystematics]
        avail = list(variProxy.calc.available())
        enable = [ vari for vari in avail if vari != "nominal" and any(vari.startswith(ena) for ena in enableSystematics) ]
        if enable:
            for opWithSyst in variProxy.brMapMap[variProxy.withSystName].values():
                opWithSyst.variations = enable ## fine, just (re)creataed by _initFromCalc
            logger.info("Enabled systematic variations for MET: {0}".format(" ".join(enable)))

def forceDefine(arg, selection):
    """ Force the definition of an expression as a column at a selection stage

    Use only for really computation-intensive operations that need to be precalculated

    :param arg: expression to define as a column
    :param selection: :py:class:`~bamboo.plots.Selection` for which the expression should be defined
    """
    if arg is None:
        raise RuntimeError("Trying to define None. If a correction calculator product was passed: has the calculator been added and configured?")
    from .treeoperations import adaptArg
    op = adaptArg(arg)
    if not op.canDefine:
        raise RuntimeError("Cannot define {0!r}".format(op))
    return selection._fbe.selDFs[selection.name].define(op)

def makePileupWeight(puWeightsFile, numTrueInteractions, variation="Nominal", systName=None, nameHint=None):
    """ Construct a pileup weight for MC, based on the weights in a JSON file

    :param puWeightsFile: path of the JSON file with weights (binned in NumTrueInteractions)
    :param numTrueInteractions: expression to get the number of true interactions (Poissonian expectation value for an event)
    :param systName: name of the associated systematic nuisance parameter
    """
    from . import treefunctions as op
    from .treeoperations import ScaleFactorWithSystOp
    paramVType = "Parameters::value_type::value_type"
    puArgs = op.construct("Parameters", (op.initList("std::initializer_list<{0}>".format(paramVType), paramVType,
        (op.initList(paramVType, "float", (op.extVar("int", "BinningVariable::NumTrueInteractions"), numTrueInteractions)),)),))
    puWFun = op.define("ILeptonScaleFactor", 'const ScaleFactor <<name>>{{"{0}"}};'.format(puWeightsFile), nameHint=nameHint)
    expr = puWFun.get(puArgs, op.extVar("int", variation))
    if systName and variation == "Nominal": ## wrap
        expr._parent = ScaleFactorWithSystOp(expr._parent, systName)
    return expr

def makeMultiPrimaryDatasetTriggerSelection(sampleName, datasetsAndTriggers):
    """ Construct a selection that prevents processing multiple times (from different primary datasets)

    If an event is passes triggers for different primary datasets, it will be taken
    from the first of those (i.e. the selection will be 'passes one of the triggers that
    select it for this primary dataset, and not for any of those that come before in the
    input dictionary).

    :param sampleName: sample name
    :param datasetsAndTriggers: a dictionary ``{primary-dataset, set-of-triggers}``, where
        the key is either a callable that takes a sample name and returns true in case
        it originates from the corresponding primary datasets, or a string that is
        the first part of the sample name in that case. The value (second item) can be
        a single expression (e.g. a trigger flag, or an OR of them), or a list of those
        (in which case an OR-expression is constructed from them).
    :returns: an expression to filter the events in the sample with given name

    :Example:

    >>> if not self.isMC(sample):
    >>>     trigSel = noSel.refine("trigAndPrimaryDataset",
    >>>         cut=makeMultiPrimaryDatasetTriggerSelection(sample, {
    >>>               "DoubleMuon" : [ t.HLT.Mu17_TrkIsoVVL_Mu8_TrkIsoVVL, t.HLT.Mu17_TrkIsoVVL_TkMu8_TrkIsoVVL ],
    >>>               "DoubleEG"   : t.HLT.Ele23_Ele12_CaloIdL_TrackIdL_IsoVL_DZ,
    >>>               "MuonEG"     : [ t.HLT.Mu23_TrkIsoVVL_Ele12_CaloIdL_TrackIdL_IsoVL, t.HLT.Mu8_TrkIsoVVL_Ele23_CaloIdL_TrackIdL_IsoVL ]
    >>>               }))
    """
    # python3.6+ dictionaries keep insertion order (implementation detail in cpython 3.6, part of the language spec since 3.7)
    # see https://docs.python.org/3/library/stdtypes.html#mapping-types-dict
    from . import treefunctions as op
    inDSAndTrigSel = [
            ((sampleName.startswith(dsArg) if isinstance(dsArg, str) else dsArg(sampleName)),
             (op.OR(*trigArg) if hasattr(trigArg, "__iter__") else trigArg))
        for dsArg, trigArg in ( datasetsAndTriggers.items() if isinstance(datasetsAndTriggers, dict) else datasetsAndTriggers) ]
    sels_not = []
    trigSel = None
    for inThisDataset, dsTrigSel in inDSAndTrigSel:
        if not inThisDataset:
            sels_not.append(dsTrigSel)
        else:
            trigSel = dsTrigSel
            break
    if trigSel is None:
        raise RuntimeError("Sample name {0} matched none of the primary datasets".format(sampleName))
    if len(sels_not) == 0: ## first
        return trigSel
    else:
        return op.AND(op.NOT(op.OR(*sels_not)), trigSel)

def configureRochesterCorrection(variProxy, paramsFile, isMC=False, backend=None, uName=""):
    """ Apply the Rochester correction for muons

    :param variProxy: muon variatons proxy, e.g. ``tree.._Muon`` for NanoAOD
    :param paramsFile: path of the text file with correction parameters
    :param isMC: MC or not
    :param backend: backend pointer (returned from :py:meth:`~bamboo.analysismodules.HistogramsModule.prepareTree`)
    :param uName: unique name for the correction calculator (sample name is a safe choice)
    """
    from bamboo.treefunctions import _to, _tp
    aMu = variProxy.orig[0]
    args = [ getattr(aMu, comp).op.arg for comp in ("pt", "eta", "phi", "mass", "charge", "nTrackerLayers") ]
    if isMC:
        args += [ aMu.genPart._idx.arg, variProxy._parent.GenPart[0].pt.op.arg ]
        evt = variProxy._parent
        args.append((evt.run<<20) + (evt.luminosityBlock<<10) + evt.event + 169)
    else:
        args += [ _to.ExtVar("ROOT::VecOps::RVec<Int_t>", "ROOT::VecOps::RVec<Int_t>{}"),
                  _to.ExtVar("ROOT::VecOps::RVec<float>", "ROOT::VecOps::RVec<float>{}"),
                  _tp.makeConst(0, "unsigned") ] ## no seed
    ## load necessary library and header(s)
    from . import treefunctions as op # first to load default headers/libs, if still needed
    from .root import loadRochesterCorrectionCalculator, gbl
    loadRochesterCorrectionCalculator()
    ## define calculator and initialize
    roccorName = backend.symbol("RochesterCorrectionCalculator <<name>>{{}}; // for {0}".format(uName), nameHint="bamboo_roccorCalc{0}".format("".join(c for c in uName if c.isalnum())))
    variProxy._initCalc(op.extVar("RochesterCorrectionCalculator", roccorName), calcHandle=getattr(gbl, roccorName), args=args)
    if not os.path.exists(paramsFile):
        raise ValueError("File {0} not found".format(paramsFile))
    variProxy.calc.setRochesterCorrection(paramsFile)
    variProxy._initFromCalc()
    assert variProxy.calcProd
