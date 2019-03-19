import base64
import json
import logging
logger = logging.getLogger(__name__)
import os
import os.path
import requests

class JetDatabaseCache(object):
    def __init__(self, name, service="https://api.github.com/repos", repository=None, branch="master", cachedir=None):
        if cachedir is None:
            from .analysisutils import bamboo_cachedir
            cachedir = bamboo_cachedir
        self.cachedir = os.path.join(cachedir, name)
        if not os.path.isdir(self.cachedir):
            os.makedirs(self.cachedir)
        self.statusFile = os.path.join(self.cachedir, "status.json")
        self.service = service
        self.repository = repository
        self.branch = branch
        self._status = None
        self._baseUrl = None

    def _get_master_sha(self):
        headers = {}
        if "branches_etag" in self._status and "sha" in self._status:
            headers["If-None-Match"] = self._status["branches_etag"]
        r_branches = requests.get("{0}/git/refs/heads".format(self._baseUrl), headers=headers)
        if r_branches.status_code == 304:
            return self._status["sha"]
        else:
            r_master = next(itm for itm in r_branches.json() if itm["ref"] == "refs/heads/master")
            self._status["branches_etag"] = r_branches.headers["ETag"]
            return r_master["object"]["sha"]

    def __enter__(self):
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

        masterSHA = self._get_master_sha()
        ## Update the *tree* of the (fetched tags of) textFiles if necessary
        if "sha" not in self._status or self._status["sha"] != masterSHA:
            logger.debug("Updating root of {0} at {1}".format(self.repository, masterSHA))
            r_rootTree = requests.get("{0}/git/trees/{1}".format(self._baseUrl, masterSHA)).json()
            r_TF = next(itm for itm in r_rootTree["tree"] if itm["path"] == "textFiles")
            if "textFiles" not in self._status:
                self._status["textFiles"] = {"tree":{}}
            statTF = self._status["textFiles"]
            if statTF.get("sha") != r_TF["sha"]:
                stTags = statTF["tree"]
                logger.debug("Updating 'textFiles' of {0} at {1} ({2})".format(self.repository, masterSHA, r_TF["sha"]))
                r_tagsTree = requests.get("{0}/git/trees/{1}".format(self._baseUrl, r_TF["sha"])).json()
                for r_tg in r_tagsTree["tree"]:
                    if r_tg["path"] not in stTags:
                        logger.debug("New tag in {0}: {1}".format(self.repository, r_tg["path"]))
                        stTags[r_tg["path"]] = {"sha": r_tg["sha"]}
                    else:
                        statTag = stTags[r_tg["path"]]
                        if r_tg["sha"] != statTag["sha"]:
                            if "tree" in statTag:
                                self._updateTag(r_tg["path"], r_tg["sha"])
                            else: ## don't trigger fetch
                                statTag["sha"] = r_tg["sha"]
                toRemove = []
                for tagName, statTag in stTags.items():
                    if not any(r_tg["path"] == tagName for r_tg in r_tagsTree["tree"]):
                        logger.debug("Tag {0} was removed from {1}".format(tagname, self.repository))
                        tagDir = os.path.join(self.cachedir, tag)
                        if os.path.isdir(tagDir):
                            import shutil
                            logger.debug("Should remove cache directory {0}".format(tagDir))
                            #shutil.rmtree(tagDir)
                        toRemove.append(tagNamE)
                for tagName in toRemove:
                    del stTags[tagName]
                statTF["sha"] = r_TF["sha"]
            self._status["sha"] = masterSHA
        else:
            logger.debug("{0} tree up to date at {1}".format(self.repository, self._status["sha"]))

        return self

    def _updateTag(self, tag, sha):
        statTag = self._status["textFiles"]["tree"][tag]
        if "tree" not in statTag:
            statTag["tree"] = {}
        tagPLs = statTag["tree"]
        logger.debug("Updating 'textFiles/{0}' of {1} to {2}".format(tag, self.repository, sha))
        r_plTree = requests.get("{0}/git/trees/{1}".format(self._baseUrl, sha)).json()
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

    def getPayload(self, tag, what, jets):
        return self._getPayload(tag, "{0}_{1}_{2}.txt".format(tag, what, jets))

    def _getPayload(self, tag, fName):
        stTags = self._status["textFiles"]["tree"]
        if tag not in stTags:
            raise ValueError("Unknown tag: {0}".format(tag))
        stTag = stTags[tag]
        if "tree" not in stTag:
            self._updateTag(tag, stTag["sha"])
        if fName not in stTag["tree"]:
            raise ValueError("Unknown file: {0}/{1}".format(tag, fName))
        stPL = stTag["tree"][fName]
        if len(stPL) == 1:
            logger.debug("Getting payload for {0}/{1} ({2})".format(tag, fName, stPL["sha"]))
            r_pl = requests.get("{0}/git/blobs/{1}".format(self._baseUrl, stPL["sha"])).json()
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
            return self._getPayload(*stPL["symlink"])
        elif "path" in stPL:
            return stPL["path"]

    def __exit__(self, type, value, traceback):
        with open(self.statusFile, "w") as sFile:
            json.dump(self._status, sFile)

if __name__ == "__main__":
    logging.basicConfig(level=logging.DEBUG)
    import IPython
    with JetDatabaseCache("JECDB", repository="cms-jet/JECDatabase", cachedir=".") as jecDBCache:
        IPython.embed()
