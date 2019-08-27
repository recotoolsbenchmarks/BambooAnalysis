"""
Batch cluster tools (common part for HTCondor and slurm, potentially others)
"""

from itertools import chain
import logging
logger = logging.getLogger(__name__)
import os, os.path
import subprocess
from datetime import datetime
import time

class CommandListJob(object):
    """ Interface/base for 'backend' classes to create a job cluster/array from a list of commands (each becoming a subjob) """
    def __init__(self, commandList, workDir=None, workdir_default_pattern="batch_work"):
        self.commandList = commandList
        self.workDir = self.init_dirtowrite(workDir, default_pattern=workdir_default_pattern)
        self.workDirs = self.setupBatchDirs(self.workDir)

    ## interface methods
    def submit(self):
        """ Submit the job(s) """
        pass
    def commandOutFiles(self, command):
        """ Output files for a given command (when finished) """
        pass
    ## for monitoring
    @property
    def status(self):
        """ overall job status summary """
        pass
    def statuses(self):
        """ list of subjob statuses (numeric, using indices in CondorJobStatus) """
        pass
    def updateStatuses(self):
        """ Update the subjob status cache (if used in the implementation) """
        pass
    def subjobStatus(self, i):
        """ get status of subjob """
        pass
    def commandStatus(self, command):
        """ get the status of the jobs corresponding to one of the commands """
        pass
    def getLogFile(self, command):
        """ get path to log file corresponding to the given command """

    ## helper methods
    @staticmethod
    def init_dirtowrite(given=None, default_pattern="batch"):
        """ Initialisation helper: check that the directory does not exist, then create it

        If None is given, a default ("$(pwd)/{default_pattern}", "$(pwd)/default_pattern_0" etc.) is used
        """
        if given is None:
            test_dir = os.path.join(os.getcwd(), default_pattern)
            if os.path.exists(test_dir):
                i = 0
                while os.path.exists(test_dir) and i < 10:
                    test_dir = os.path.join(os.getcwd(), "{0}_{1:d}".format(default_pattern, i))
                    i += 1
        else:
            test_dir = given
        if os.path.exists(test_dir):
            raise RuntimeError("Directory {0} exists".format(test_dir))

        os.makedirs(test_dir)

        return os.path.abspath(test_dir) ## make sure we keep absolute paths

    @staticmethod
    def setupBatchDirs(workDir):
        """ Create up the working directories (input, output, logs) under workDir """
        dirs = {
              "in"  : os.path.join(workDir, "input")
            , "out" : os.path.join(workDir, "output")
            , "log" : os.path.join(workDir, "logs")
            }
        for subdir in dirs.values():
            os.mkdir(subdir)
        logger.info("Created working directories under {0}".format(workDir))
        return dirs

class SplitAggregationTask(object):
    """ Task split in commands, whose results can be merged """
    def __init__(self, commandList, finalizeAction=None):
        self._jobCluster = None
        self.commandList = commandList
        self.finalizeAction = finalizeAction
        self.failedCommands = set()
    @property
    def jobCluster(self):
        return self._jobCluster
    @jobCluster.setter
    def jobCluster(self, jobCluster):
        if self._jobCluster:
            raise Exception("SplitAggregationTask.jobCluster has already been set")
        self._jobCluster = jobCluster
        if self.finalizeAction: ## pass it on
            self.finalizeAction.jobCluster = jobCluster

    def tryFinalize(self, failedStatuses):
        """Return True if all jobs have finished running (either failed or completed)"""
        if self.jobCluster:
            self.failedCommands.update(cmd for cmd in self.commandList if self.jobCluster.commandStatus(cmd).upper() in failedStatuses)
            if len(self.failedCommands) == len(self.commandList):
                return True
            if all(self.jobCluster.commandStatus(cmd).upper() == "COMPLETED" for cmd in self.commandList):
                if self.finalizeAction:
                    self.finalizeAction.perform()
                return True
        return False

class Action(object):
    """ interface for job finalization """
    def perform(self, dryRun=False):
        """ interface method """
        return False

class HaddAction(Action):
    """ merge results with hadd
    """
    def __init__(self, commandList, outDir, options=None):
        self.jobCluster = None
        self.commandList = commandList
        self.outDir = outDir
        self.options = options if options is not None else []
        super(HaddAction, self).__init__()
    def perform(self, dryRun=False):
        """ Perform actual hadd/move. If dryRun is True, do not call the commands but just print them."""
        if len(self.commandList) == 0:   ## nothing
            return
        elif len(self.commandList) == 1: ## move
            cmd = self.commandList[0]
            for outf in self.jobCluster.commandOutFiles(cmd):
                action = ["mv", outf, self.outDir]
                if not dryRun:
                    logger.info("Finalization: moving output file {0} to {1}".format(outf, self.outDir))
                    try:
                        subprocess.check_call(action)
                    except subprocess.CalledProcessError:
                        logger.error("Failed to run the finalization command:")
                        logger.error(" ".join(action))
                else:
                    logger.error("Not performing finalization with command:")
                    logger.error(" ".join(action))
        else:                            ## merge
            ## collect for each output file name which jobs produced one
            filesToMerge = dict()
            for cmd in self.commandList:
                for outf in self.jobCluster.commandOutFiles(cmd):
                    outfb = os.path.basename(outf)
                    if outfb not in filesToMerge:
                        filesToMerge[outfb] = []
                    filesToMerge[outfb].append(outf)

            for outfb, outfin in filesToMerge.items():
                fullout = os.path.join(self.outDir, outfb)
                action = ["hadd"]+self.options+[fullout]+outfin
                if not dryRun:
                    logger.info("Finalization: merging {0} to {1}".format(", ".join(outfin), fullout))
                    try:
                        subprocess.check_call(action)
                    except subprocess.CalledProcessError:
                        logger.error("Failed to run the finalization command:")
                        logger.error(" ".join(action))
                else:
                    logger.error("Not performing finalization with command:")
                    logger.error(" ".join(action))

class TasksMonitor(object):
    """ Monitor a number of tasks and the associated job clusters """
    def __init__(self, jobs=[], tasks=[], interval=120, allStatuses=None, activeStatuses=None, failedStatuses=None, completedStatus=None):
        """ Constructor

        jobs: job clusters to monitor
        tasks: tasks to monitor (and finalize)
        interval: number of seconds between status checks
        """
        self.jobs = set()
        self.tasks = set()
        self.activeTasks = set()
        self.add(jobs, tasks)
        self.interval = interval
        self.allStatuses = list(allStatuses) if allStatuses else []
        self.activeStatuses = list(activeStatuses) if activeStatuses else []
        self.failedStatuses = list(failedStatuses) if failedStatuses else []
        self.completedStatus = completedStatus
    def add(self, jobs, tasks):
        """ add jobs and tasks """
        self.jobs.update(jobs)
        self.tasks.update(tasks)
        assert sum(len(j.commandList) for j in self.jobs) == sum(len(t.commandList) for t in self.tasks) ## assumption: same set of commands
        self.activeTasks.update(tasks)
    @staticmethod
    def makeStats(statuses, allStatuses):
        histo = [ 0 for st in allStatuses ]
        for jSt in statuses:
            histo[jSt] += 1
        return histo
    @staticmethod
    def formatStats(stats, allStatuses):
        return "{0} / {1:d} Total".format(", ".join("{0:d} {1}".format(n,nm) for n,nm in zip(stats, allStatuses)), sum(stats))

    def _tryFinalize(self):
        """ try to finalize tasks """
        finalized = []
        for t in self.activeTasks:
            if t.tryFinalize([ self.allStatuses[i] for i in self.failedStatuses ]):
                finalized.append(t) ## don't change while iterating
        for ft in finalized:
            self.activeTasks.remove(ft)
        if len(finalized) > 0:
            logger.info("Finalized {0:d} tasks, {1:d} remaining".format(len(finalized), len(self.activeTasks)))

    def collect(self, wait=None):
        """ Wait for the jobs to finish, then finalize tasks.
        Returns True if every job succeeded, otherwise return False.
        """
        shouldDoPostProc = True
        for j in self.jobs:
            j.updateStatuses()
        self._tryFinalize()
        if len(self.activeTasks) > 0:
            logger.info("Waiting for jobs to finish...")
            nJobs = sum( len(j.commandList) for j in self.jobs )
            if wait:
                time.sleep(wait)
            stats = self.makeStats(chain.from_iterable(j.statuses() for j in self.jobs), self.allStatuses)
            while len(self.activeTasks) > 0 and sum(stats[sa] for sa in self.activeStatuses) > 0:
                time.sleep(self.interval)
                prevStats = stats
                stats = self.makeStats(chain.from_iterable(j.statuses() for j in self.jobs), self.allStatuses)
                logger.info("[ {0} :: {1} ]".format(datetime.now().strftime("%H:%M:%S"), self.formatStats(stats, self.allStatuses)))
                if stats[self.completedStatus] > prevStats[self.completedStatus]:
                    self._tryFinalize()
        # Check for failed tasks, print manual recovering actions:
        for task in self.tasks:
            if task.failedCommands:
                shouldDoPostProc = False
                logger.error("The following commands failed to run:")
                for cmd in task.failedCommands:
                    logger.error(cmd)
                logger.error("The job logs are in:")
                for cmd in task.failedCommands:
                    logger.error(task.jobCluster.getLogFile(cmd))
                if task.finalizeAction:
                    logger.error("The corresponding finalization command is:")
                    task.finalizeAction.perform(dryRun=True)
        return shouldDoPostProc

def splitInChunks(theList, nChunks=None, chunkLength=None):
    if ( nChunks is None ) == ( chunkLength is None ):
        raise ValueError("One and only one of nChunks or chunkLength should be specified, got {0} and {1}".format(nChunks, chunkLength))
    from itertools import islice
    N = len(theList)
    if nChunks is not None:
        import math
        chunkLength = int(math.ceil(1.*N/nChunks))
    for iStart, iStop in zip(range(0, len(theList), chunkLength), range(chunkLength, len(theList)+chunkLength, chunkLength)):
        yield islice(theList, iStart, min(iStop,N))

def writeFileList(fileList, outName):
    if os.path.isfile(outName):
        with open(outName) as ffile:
            ffiles = [ ln.strip() for ln in ffile if len(ln.strip()) > 0 ]
        if sorted(ffiles) == sorted(fileList):
            return
        else:
            logger.warning("Overwriting {0} with a new list of {1:d} files".format(outName, len(fileList)))
    with open(outName, "w") as nfile:
        nfile.write("\n".join(fileList))
