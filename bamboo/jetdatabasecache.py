from contextlib import contextmanager
import base64
import json
import logging
logger = logging.getLogger(__name__)
import os
import os.path
import requests

class StatusPermissionError(Exception):
    """ Exception for trying to write to locked status file (for internal use) """
    def __init__(self, statusFile, *args, **kwargs):
        self.statusFile = statusFile
        super(StatusPermissionError, self).__init__(*args, **kwargs)
    def __str__(self):
        return "No permission to update status file {0}".format(self.statusFile)

@contextmanager
def sessionWithResponseChecks():
    def response_check(resp, *args, **kwargs):
        try:
            resp.raise_for_status()
        except requests.exceptions.HTTPError as ex:
            logger.exception(f"Problem with request '{resp.request.method} {resp.request.url}'")
    session = requests.Session()
    session.hooks = {
        "response" : response_check
        }
    yield session

class JetDatabaseCache(object):
    def __init__(self, name, service="https://api.github.com/repos", repository=None, branch="master", cachedir=None, mayWrite=True, session=None):
        self._mayWrite = mayWrite
        if cachedir is None:
            from .analysisutils import bamboo_cachedir
            cachedir = bamboo_cachedir
        self.cachedir = os.path.join(cachedir, name)
        self.statusFile = os.path.join(self.cachedir, "status.json")
        self._statusLockName = "{0}.lock".format(self.statusFile)
        self.service = service
        self.repository = repository
        self.branch = branch
        self._status = None
        self._baseUrl = None
        if session:
            self._init(session=session)
        else:
            with sessionWithResponseChecks() as session:
                self._init(session=session)

    @contextmanager
    def _statusLockAndSave(self, expires=None):
        if not self._mayWrite:
            raise StatusPermissionError(self.statusFile)
        import time
        from datetime import datetime
        while True:
            if not os.path.exists(self._statusLockName):
                if not os.path.isdir(self.cachedir):
                    os.makedirs(self.cachedir, exist_ok=True)
                mypid = "{0}:{1:d}\n".format(os.uname()[1], os.getpid())
                try:
                    with open(self._statusLockName, "x") as lf:
                        lf.write(mypid)
                except OSError as ex:
                    logger.debug("Failed to acquire lock ({0}), waiting a bit to try again".format(str(ex)))
                except FileExistsError as ex:
                    logger.debug("Failed to acquire lock ({0}), waiting a bit to try again".format(str(ex)))
                else:
                    logger.debug("Acquired write lock for {0}".format(self.statusFile))
                    yield
                    with open(self.statusFile, "w") as sFile:
                        json.dump(self._status, sFile)
                    os.remove(self._statusLockName)
                    logger.debug("Released write lock for {0}".format(self.statusFile))
                    break
                time.sleep(5)
            elif expires is not None and ( ( datetime.now() - datetime.fromtimestamp(os.path.getctime(self._statusLockName)) ) > expires ):
                logger.warning("Lock file expired, removing")
                os.remove(self._statusLockName)
            else:
                logger.debug("{0} locked, waiting another 5s".format(self.statusFile))
                time.sleep(5) # wait my turn

    def _get_master_sha(self, session=None):
        headers = {}
        if "branches_etag" in self._status and "sha" in self._status:
            headers["If-None-Match"] = self._status["branches_etag"]
        r_branches = session.get("{0}/git/refs/heads".format(self._baseUrl), headers=headers)
        if r_branches.status_code == 304:
            return self._status["sha"]
        else:
            r_master = next(itm for itm in r_branches.json() if itm["ref"] == "refs/heads/master")
            self._status["branches_etag"] = r_branches.headers["ETag"]
            return r_master["object"]["sha"]

    def _init(self, session=None):
        if os.path.exists(self.statusFile):
            with open(self.statusFile) as sFile:
                self._status = json.load(sFile)
            if "service" not in self._status:
                self._status["service"] = self.service
            elif self._status["service"] != self.service:
                raise RuntimeError("Service {0} doesn't match with status file value {1}".format(self.service, self._status["service"]))
            if "repository" not in self._status:
                self._status["repository"] = self.repository
            elif self._status["repository"] != self.repository:
                raise RuntimeError("Repository name {0} doesn't match with status file value {1}".format(self.repository, self._status["repository"]))
            self._baseUrl = "/".join((self.service.rstrip("/"), self.repository.lstrip("/"))).rstrip("/")
        else:
            self._status = {
                "service": self.service,
                "repository": self.repository,
                "branch": self.branch
                }
            self._baseUrl = "/".join((self.service.rstrip("/"), self.repository.lstrip("/"))).rstrip("/")

        masterSHA = self._get_master_sha(session=session)
        ## Update the *tree* of the (fetched tags of) textFiles if necessary
        if "sha" not in self._status or self._status["sha"] != masterSHA:
            logger.debug("Updating root of {0} at {1}".format(self.repository, masterSHA))
            with self._statusLockAndSave():
                r_rootTree = session.get("{0}/git/trees/{1}".format(self._baseUrl, masterSHA)).json()
                r_TF = next(itm for itm in r_rootTree["tree"] if itm["path"] == "textFiles")
                if "textFiles" not in self._status:
                    self._status["textFiles"] = {"tree":{}}
                statTF = self._status["textFiles"]
                if statTF.get("sha") != r_TF["sha"]:
                    stTags = statTF["tree"]
                    logger.debug("Updating 'textFiles' of {0} at {1} ({2})".format(self.repository, masterSHA, r_TF["sha"]))
                    r_tagsTree = session.get("{0}/git/trees/{1}".format(self._baseUrl, r_TF["sha"])).json()
                    for r_tg in r_tagsTree["tree"]:
                        if r_tg["path"] not in stTags:
                            logger.debug("New tag in {0}: {1}".format(self.repository, r_tg["path"]))
                            stTags[r_tg["path"]] = {"sha": r_tg["sha"]}
                        else:
                            statTag = stTags[r_tg["path"]]
                            if r_tg["sha"] != statTag["sha"]:
                                if "tree" in statTag:
                                    self._updateTag(r_tg["path"], r_tg["sha"], session=session)
                                else: ## don't trigger fetch
                                    statTag["sha"] = r_tg["sha"]
                    toRemove = []
                    for tagName, statTag in stTags.items():
                        if not any(r_tg["path"] == tagName for r_tg in r_tagsTree["tree"]):
                            logger.debug("Tag {0} was removed from {1}".format(tagName, self.repository))
                            tagDir = os.path.join(self.cachedir, tagName)
                            if os.path.isdir(tagDir):
                                import shutil
                                logger.debug("Should remove cache directory {0}".format(tagDir))
                                #shutil.rmtree(tagDir)
                            toRemove.append(tagName)
                    for tagName in toRemove:
                        del stTags[tagName]
                    statTF["sha"] = r_TF["sha"]
                self._status["sha"] = masterSHA
        else:
            logger.debug("{0} tree up to date at {1}".format(self.repository, self._status["sha"]))

        return self

    def _updateTag(self, tag, sha, session=None):
        statTag = self._status["textFiles"]["tree"][tag]
        if "tree" not in statTag:
            statTag["tree"] = {}
        tagPLs = statTag["tree"]
        logger.debug("Updating 'textFiles/{0}' of {1} to {2}".format(tag, self.repository, sha))
        r_plTree = session.get("{0}/git/trees/{1}".format(self._baseUrl, sha)).json()
        for r_pl in r_plTree["tree"]:
            if r_pl["path"] not in tagPLs:
                logger.debug("New payload in {0}/{1}: {2} ({3})".format(self.repository, tag, r_pl["path"], r_pl["sha"]))
                tagPLs[r_pl["path"]] = {"sha": r_pl["sha"]}
            else:
                statPL = tagPLs[r_pl["path"]]
                if r_pl["sha"] != statPL["sha"]:
                    logger.debug("Updated payload in {0}/{1}: {2} ({3})".format(self.repository, tag, r_pl["path"], r_pl["sha"]))
                    if "path" in statPL and os.path.isfile(statPL["path"]):
                        logger.debug("Removing outdated cached payload {0}".format(statPL["path"]))
                        os.remove(statPL["path"])
                    tagPLs[r_pl["path"]] = {"sha": r_pl["sha"]}
        toRemove = []
        for plName, statPL in tagPLs.items():
            if not any(r_pl["path"] == plName for r_pl in r_plTree["tree"]):
                logger.debug("Payload {0}/{1} was removed upstream".format(tag, plName))
                if "path" in statPL and os.path.isfile(statPL["path"]):
                    logger.debug("Should remove outdated cached payload {0}".format(statPL["path"]))
                    #os.remove(statPL["path"])
                toRemove.append(plName)
        for plName in toRemove:
            del tagPLs[plName]
        statTag["sha"] = sha

    def getPayload(self, tag, what, jets, session=None):
        try:
            plName = "{0}_{1}_{2}.txt".format(tag, what, jets)
            if not session:
                with sessionWithResponseChecks() as session:
                    return self._getPayload(tag, plName, session=session)
            else:
                return self._getPayload(tag, plName, session=session)
        except StatusPermissionError as ex:
            name = os.path.basename(self.cachedir)
            cachedir = os.path.dirname(self.cachedir)
            vName = name.lower()
            cmds = ['import logging', 'logging.basicConfig(level=logging.DEBUG)',
                    'from bamboo.jetdatabasecache import JetDatabaseCache',
                    '{0} = JetDatabaseCache("{1}", service="{2}", repository="{3}", branch="{4}", cachedir="{5}")'.format(vName, name, self.service, self.repository, self.branch, cachedir),
                    '{0}.getPayload("{1}", "{2}", "{3}")'.format(vName, tag, what, jets)
                   ]
            raise RuntimeError("\n>>> ".join(["Failed to update {0} (no permission).\nPlease update the cache by running in non-distributed mode (e.g. with --maxFile=1), or manually with".format(ex.statusFile)]+cmds))

    def _getPayload(self, tag, fName, session=None):
        stTags = self._status["textFiles"]["tree"]
        if tag not in stTags:
            raise ValueError("Unknown tag: {0}".format(tag))
        stTag = stTags[tag]
        if "tree" not in stTag:
            with self._statusLockAndSave():
                self._updateTag(tag, stTag["sha"], session=session)
        if fName not in stTag["tree"]:
            raise ValueError("Unknown file: {0}/{1}".format(tag, fName))
        stPL = stTag["tree"][fName]
        if len(stPL) == 1:
            logger.debug("Getting payload for {0}/{1} ({2})".format(tag, fName, stPL["sha"]))
            with self._statusLockAndSave():
                r_pl = session.get("{0}/git/blobs/{1}".format(self._baseUrl, stPL["sha"])).json()
                if r_pl["encoding"] == "base64":
                    content = base64.b64decode(r_pl["content"])
                else:
                    raise RuntimeError("Unknown encoding: {0}".format(r_pl["encoding"])) 
                if content.startswith(b"../") and len(content.split(b"\n")) == 1:
                    _,actTag,actFName = content.decode().strip().split(os.sep)
                    stPL["symlink"] = (actTag, actFName)
                else:
                    tagDir = os.path.join(self.cachedir, tag)
                    if not os.path.exists(tagDir):
                        os.makedirs(tagDir)
                    plFName = os.path.join(tagDir, fName)
                    with open(plFName, "wb") as plF:
                        plF.write(content)
                    stPL["path"] = plFName
                    logger.debug("Saved as {0}".format(plFName))
        if "symlink" in stPL:
            logger.debug("Payload {0}/{1} is symlinked to {2}/{3}".format(tag, fName, *stPL["symlink"]))
            return self._getPayload(*stPL["symlink"], session=session)
        elif "path" in stPL:
            return stPL["path"]

def checkCMS_CLI():
    """ Command-line script to update and check the status of the JEC/JER database caches """
    import argparse
    parser = argparse.ArgumentParser(description="Update or explore the JEC/JER database caches")
    parser.add_argument("-q", "--quiet", action="store_true", help="Run in quiet mode")
    parser.add_argument("-i", "--interactive", action="store_true", help="Open an IPython shell to explore the database caches")
    parser.add_argument("--readonly", action="store_true", help="Open the caches in readonly mode (they will not be able to automatically update, in this case)")
    parser.add_argument("--cachedir", type=str, help="Alternative cache root directory")
    args = parser.parse_args()

    logging.basicConfig(level=(logging.INFO if args.quiet else logging.DEBUG))
    mayWrite = not args.readonly
    with sessionWithResponseChecks() as session:
        logger.info("Loading JECDatabase from 'cms-jet/JECDatabase', with {0}".format(f"cache dir in {args.cachedir}" if args.cachedir else "default cache dir"))
        jecDBCache = JetDatabaseCache("JECDatabase", repository="cms-jet/JECDatabase", cachedir=args.cachedir, mayWrite=(not args.readonly))
        logger.info("Loading JRDatabase from 'cms-jet/JRDatabase', with {0}".format(f"cache dir in {args.cachedir}" if args.cachedir else "default cache dir"))
        jrDBCache = JetDatabaseCache("JRDatabase", repository="cms-jet/JRDatabase", cachedir=args.cachedir, mayWrite=(not args.readonly))
        if args.interactive:
            import IPython
            IPython.embed(header='JEC and JR database caches are available in `jecDBCache` and `jrDBCache`\nExample: pl = jecDBCache.getPayload("Summer19UL17_RunE_V1_SimpleL1_DATA", "L1FastJet", "AK4PFchs")')
