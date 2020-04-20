"""
Slurm tools (based on previous condorhelpers and cp3-llbb/CommonTools condorSubmitter and slurmSubmitter)
"""
__all__ = ("CommandListJob", "jobsFromTasks", "makeTasksMonitor")

from itertools import chain
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

SlurmJobStatus = ["PENDING", "RUNNING", "COMPLETED", "FAILED", "COMPLETING", "CONFIGURING", "CANCELLED", "BOOT_FAIL", "NODE_FAIL", "PREEMPTED", "RESIZING", "SUSPENDED", "TIMEOUT", "OUT_OF_MEMORY", "unknown"]
SlurmJobStatus_active = ["CONFIGURING", "COMPLETING", "PENDING", "RUNNING", "RESIZING", "SUSPENDED"]
SlurmJobStatus_failed = ["FAILED", "TIMEOUT", "CANCELLED", "OUT_OF_MEMORY"]
SlurmJobStatus_completed = "COMPLETED"

class CommandListJob(CommandListJobBase):
    """
    Helper class to create a slurm job array from a list of commands (each becoming a task in the array)
    
    Default work directory will be $(pwd)/slurm_work, default output pattern is "*.root"
    """
    default_cfg_opts = {
          "batchScriptsFilename" : "slurmSubmission.sh"
        , "stageoutFiles"        : ["*.root"]
        , "inputParamsNames"     : ["taskcmd"]
        , "useJobArray"          : True
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

    def _arrayIndex(self, command):
        return self.commandList.index(command)+1

    def _commandOutDir(self, command):
        """ Output directory for a given command """
        return os.path.join(self.workDirs["out"], str(self._arrayIndex(command)))

    def commandOutFiles(self, command):
        """ Output files for a given command """
        import fnmatch
        cmdOutDir = self._commandOutDir(command)
        if not os.path.isdir(cmdOutDir):
            os.mkdir(cmdOutDir)
        return list( os.path.join(cmdOutDir, fn) for fn in os.listdir(cmdOutDir)
                if any( fnmatch.fnmatch(fn, pat) for pat in self.cfg.stageoutFiles) )

    def submit(self):
        """ Submit the job to slurm """
        sbatchOpts = [
              "--partition={}".format(self.cfg.sbatch_partition)
            , "--qos={}".format(self.cfg.sbatch_qos)
            ]+(list(self.cfg.sbatch_additionalOptions) if hasattr(self.cfg, "sbatch_additionalOptions") and self.cfg.sbatch_additionalOptions else [])
        logger.info("Submitting an array of {0:d} jobs to slurm".format(len(self.commandList)))
        logger.debug("sbatch {0} {1}".format(" ".join(sbatchOpts), self.slurmScript))
        result = subprocess.check_output(["sbatch"]+sbatchOpts+[self.slurmScript]).decode()

        self.clusterId = next(tok for tok in reversed(result.split()) if tok.isdigit())
        
        ## save to file in case
        with open(os.path.join(self.workDirs["in"], "cluster_id"), "w") as f:
            f.write(self.clusterId)

        logger.info("Submitted, job ID is {}".format(self.clusterId))

    def cancel(self):
        """ Cancel all the jobs using scancel """

        subprocess.check_call(["scancel", self.clusterId])

    def statuses(self, update=True):
        """ list of subjob statuses (numeric, using indices in SlurmJobStatus) """
        if update:
            try:
                self.updateStatuses()
            except Exception as ex:
                logger.error("Exception while updating statuses (will reuse previous): {0!s}".format(ex))
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
                try:
                    sacctCmdArgs = ["sacct", "-X", "-n", "--format", "State%20", "-j", subjobId]
                    ret = subprocess.check_output(sacctCmdArgs).decode().strip()
                    if "\n" in ret:
                        raise AssertionError("More than one line in sacct... there's something wrong")
                    if len(ret) != 0:
                        if "CANCELLED+" in ret:
                            # Can happen if scancel command did not have time to propagate
                            status = "CANCELLED"
                        else:
                            status = ret
                    else:
                        squeueCmdArgs = ["squeue", "-h", "-O", "state", "-j", subjobId]
                        status = subprocess.check_output(squeueCmdArgs).decode().strip()
                    if status not in SlurmJobStatus_active:
                        self._finishedTasks[subjobId] = status
                except subprocess.CalledProcessError as ex:
                    logger.warning(f"Could not update status of job {subjobId}: {ex!s}")
                    # fall back to previous status (probably PENDING or RUNNING) and try again later
                    status = self._statuses[i]
            self._statuses[i] = status

    def commandStatus(self, command):
        return self.subjobStatus(self._arrayIndex(command))

    def getLogFile(self, command):
        return os.path.join(self.workDirs["log"], "slurm-{}_{}.out".format(self.clusterId, self._arrayIndex(command)))

    def getResubmitCommand(self, failedCommands):
        sbatchOpts = [
              "--array={}".format(",".join(str(self._arrayIndex(cmd)) for cmd in failedCommands))
            , "--partition={}".format(self.cfg.sbatch_partition)
            , "--qos={}".format(self.cfg.sbatch_qos)
            ]+(list(self.cfg.sbatch_additionalOptions) if hasattr(self.cfg, "sbatch_additionalOptions") and self.cfg.sbatch_additionalOptions else [])
        return ["sbatch"]+sbatchOpts+[self.slurmScript]

    def getRuntime(self, command):
        sacctCmdArgs = ["sacct", "-n", "--format", "Elapsed", "-j", "{0}_{1:d}".format(self.clusterId, self._arrayIndex(command))]
        elapsed = subprocess.check_output(sacctCmdArgs).decode().split("\n")[0].strip()
        dys, ms = 0, 0
        if "-" not in elapsed:
            elapsed = "0-{0}".format(elapsed)
            dys, elapsed = (tk.strip() for tk in elapsed.split("-"))
            dys = int(dys)
        if "." in elapsed:
            elapsed, ms = (tk.strip() for tk in elapsed.split("."))
            ms = int(ms)
        import datetime
        return (datetime.datetime.strptime(elapsed, "%H:%M:%S")-datetime.datetime(1900, 1, 1)) + datetime.timedelta(days=dys, milliseconds=ms)

def jobsFromTasks(taskList, workdir=None, batchConfig=None, configOpts=None):
    if batchConfig:
        if not configOpts:
            configOpts = dict()
        else:
            configOpts = dict(configOpts)
        bc_c = dict(batchConfig)
        for k,cov in configOpts.items():
            if k.lower() in bc_c:
                if isinstance(cov, list):
                    cov += bc_c[k.lower()].split(", ")
                    del bc_c[k.lower()]
                else: # ini file takes preference
                    configOpts[k] = bc_c[k.lower()]
                    del bc_c[k.lower()]
        configOpts.update(bc_c)
    slurmJob = CommandListJob(list(chain.from_iterable(task.commandList for task in taskList)), workDir=workdir, configOpts=configOpts)
    for task in taskList:
        task.jobCluster = slurmJob
    return [ slurmJob ]

def makeTasksMonitor(jobs=[], tasks=[], interval=120):
    """ make a TasksMonitor for slurm jobs """
    return TasksMonitor(jobs=jobs, tasks=tasks, interval=interval
            , allStatuses=SlurmJobStatus
            , activeStatuses=[SlurmJobStatus.index(stNm) for stNm in SlurmJobStatus_active]
            , failedStatuses=[SlurmJobStatus.index(stNm) for stNm in SlurmJobStatus_failed]
            , completedStatus=SlurmJobStatus.index(SlurmJobStatus_completed)
            )
