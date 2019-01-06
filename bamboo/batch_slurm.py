"""
Slurm tools (based on previous condorhelpers and cp3-llbb/CommonTools condorSubmitter and slurmSubmitter)
"""
__all__ = ("CommandListJob", "jobsFromTasks", "makeTasksMonitor")

from itertools import chain
from contextlib import contextmanager
import logging
logger = logging.getLogger(__name__)
import os.path
import subprocess

from .batch import CommandListJob as CommandListJobBase
from .batch import TasksMonitor

try:
    from CP3SlurmUtils.Configuration import Configuration as CP3SlurmConfiguration
    from CP3SlurmUtils.SubmitWorker import SubmitWorker as slurmSubmitWorker
except ImportError as ex:
    logger.info("Could not import Configuration and slurmSubmitWorker from CP3SlurmUtils.SubmitUtils. Please run 'module load slurm/slurm_utils'")

SlurmJobStatus = ["PENDING", "RUNNING", "COMPLETED", "FAILED", "COMPLETING", "CONFIGURING", "CANCELLED", "BOOT_FAIL", "NODE_FAIL", "PREEMPTED", "RESIZING", "SUSPENDED", "TIMEOUT", "unknown"]

class CommandListJob(CommandListJobBase):
    """
    Helper class to create a slurm job array from a list of commands (each becoming a task in the array)
    
    Default work directory will be $(pwd)/slurm_work, default output pattern is "*.root"
    """
    default_cfg_opts = {
          "batchScriptsFilename" : "slurmSubmission.sh"
        , "stageoutFiles"        : ["*.root"]
        , "inputParamsNames"     : ["taskcmd"]
        , "useJobArray"          :  True
        , "payload"              : ("${taskcmd}")
        }

    def __init__(self, commandList, workDir=None, configOpts=None):
        super(CommandListJob, self).__init__(commandList, workDir=workDir, workdir_default_pattern="slurm_work")
        ##
        self.cfg = CP3SlurmConfiguration()

        self.cfg.sbatch_workdir = self.workDir
        self.cfg.inputSandboxDir = self.workDirs["in"]
        self.cfg.batchScriptsDir = self.workDir
        self.cfg.stageoutDir = os.path.join(self.workDirs["out"], "${SLURM_ARRAY_TASK_ID}")
        self.cfg.stageoutLogsDir = self.workDirs["log"]
        self.cfg.inputParams = list([cmd] for cmd in self.commandList)
        ## apply user-specified
        cfg_opts = dict(CommandListJob.default_cfg_opts)
        if configOpts:
            cfg_opts.update(configOpts)
        for k, v in cfg_opts.items():
            setattr(self.cfg, k, v)

        self.slurmScript = os.path.join(self.cfg.batchScriptsDir, self.cfg.batchScriptsFilename)
        self.clusterId = None ## will be set by submit
        self._finishedTasks = dict()
        self._statuses = ["PENDING" for cmd in self.commandList]

        ## output
        try:
            slurm_submit = slurmSubmitWorker(self.cfg, submit=False, debug=False, quiet=False)
        except Exception as ex:
            logger.error("Problem constructing slurm submit worker: {}".format(str(ex)))
            raise ex
        else:
            from contextlib import redirect_stdout
            from io import StringIO
            workerout = StringIO()
            submitLoggerFun = logger.debug
            try:
                with redirect_stdout(workerout):
                    slurm_submit()
            except Exception as ex:
                submitLoggerFun = logger.info
                raise RuntimeError("Problem in slurm_submit: {}".format(str()))
            finally:
                submitLoggerFun("==========     BEGIN slurm_submit output     ==========")
                for ln in workerout.getvalue().split("\n"):
                    submitLoggerFun(ln)
                submitLoggerFun("==========     END   slurm_submit output     ==========")

    def _commandOutDir(self, command):
        """ Output directory for a given command """
        return os.path.join(self.workDirs["out"], str(self.commandList.index(command)+1))

    def commandOutFiles(self, command):
        """ Output files for a given command """
        import fnmatch
        cmdOutDir = self._commandOutDir(command)
        return list( os.path.join(cmdOutDir, fn) for fn in os.listdir(cmdOutDir)
                if any( fnmatch.fnmatch(fn, pat) for pat in self.cfg.stageoutFiles) )

    def submit(self):
        """ Submit the job to slurm """
        logger.info("Submitting an array of {0:d} jobs to slurm".format(len(self.commandList)))
        result = subprocess.check_output(["sbatch"
            , "--partition={}".format(self.cfg.sbatch_partition)
            , "--qos={}".format(self.cfg.sbatch_qos)
            , "--wckey=cms"
            , self.slurmScript]).decode()

        self.clusterId = next(tok for tok in reversed(result.split()) if tok.isdigit())
        
        ## save to file in case
        with open(os.path.join(self.workDirs["in"], "cluster_id"), "w") as f:
            f.write(self.clusterId)

        logger.info("Submitted, job ID is {}".format(self.clusterId))

    def statuses(self, update=True):
        """ list of subjob statuses (numeric, using indices in SlurmJobStatus) """
        if update:
            self.updateStatuses()
        return [ SlurmJobStatus.index(sjst) for sjst in self._statuses ]

    @property
    def status(self):
        if all(st == self._statuses[0] for st in self._statuses):
            return self._statuses[0]
        elif any(st == "RUNNING" for st in self._statuses):
            return "RUNNING"
        elif any(st == "CANCELLED" for st in self._statuses):
            return "CANCELLED"
        else:
            return "unknown"

    def subjobStatus(self, i):
        return self._statuses[i-1]

    def updateStatuses(self):
        for i in range(len(self.commandList)):
            subjobId = "{0}_{1:d}".format(self.clusterId, i+1)
            status = "unknown"
            if subjobId in self._finishedTasks:
                status = self._finishedTasks[subjobId]
            else:
                sacctCmdArgs = ["sacct", "-n", "--format", "State", "-j", subjobId]
                ret = subprocess.check_output(sacctCmdArgs).decode().strip()
                if "\n" in ret: ## if finished
                    if len(ret.split("\n")) == 2:
                        ret = subprocess.check_output(["sacct", "-n", "--format", "State", "-j", "{}.batch".format(subjobId)]).decode().strip()
                        self._finishedTasks[subjobId] = ret
                        status = ret
                    else:
                        raise AssertionError("More than two lines in sacct... there's something wrong")
                else:
                    if len(ret) != 0:
                        status = ret
                    else:
                        squeueCmdArgs = ["squeue", "-h", "-O", "state", "-j", subjobId]
                        ret = subprocess.check_output(squeueCmdArgs).decode().strip()
                        if len(ret) != 0:
                            status = ret
                        else: # fall back to previous status (probably PENDING or RUNNING)
                            status = self._statuses[i]
            self._statuses[i] = status

    def commandStatus(self, command):
        return self.subjobStatus(self.commandList.index(command)+1)

def jobsFromTasks(taskList, workdir=None, batchConfig=None, configOpts=None):
    slurmJob = CommandListJob(list(chain.from_iterable(task.commandList for task in taskList)), workDir=workdir, configOpts=configOpts)
    for task in taskList:
        task.jobCluster = slurmJob
    return [ slurmJob ]

def makeTasksMonitor(jobs=[], tasks=[], interval=120):
    """ make a TasksMonitor for slurm jobs """
    return TasksMonitor(jobs=jobs, tasks=tasks, interval=interval
            , allStatuses=SlurmJobStatus
            , activeStatuses=[SlurmJobStatus.index(stNm) for stNm in ("CONFIGURING", "COMPLETING", "PENDING", "RUNNING", "RESIZING", "SUSPENDED")]
            , completedStatus=SlurmJobStatus.index("COMPLETED")
            )
