"""
Collect PyROOT imports, and utilities for loading libraries and headers
"""
__all__ = ("gbl",)
## workaround for https://sft.its.cern.ch/jira/browse/ROOT-10304, which affects only ROOT 6.18/XX
import subprocess
_rootVersion = subprocess.check_output(["root-config", "--version"]).decode().strip()
_rootVersion_split = _rootVersion.split
if _rootVersion.split("/")[0].split(".")[1] == "18":
    import ROOT as gbl
else:
    from cppyy import gbl
