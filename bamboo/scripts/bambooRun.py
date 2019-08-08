#!/usr/bin/env python
"""
Command-line script: run an analysis module

Loads a class from a module (import path or file), constructs an instance with the remaining arguments to parse, and runs it
"""

import argparse
import importlib.util
import logging
import os.path

def main():
    parser = argparse.ArgumentParser(description="Run an analysis module", add_help=False)
    parser.add_argument("-m", "--module", type=str, default="bamboo.analysismodules:AnalysisModule", help="Module to run (format: modulenameOrPath[:classname])")
    parser.add_argument("--help", "-h", action="store_true", help="Print this help message and exit")
    parser.add_argument("-v", "--verbose", action="store_true", help="Run in verbose mode")
    args, remArgs = parser.parse_known_args()
    ## pass arguments on
    remArgs.append("--module={0}".format(args.module))
    if args.help:
        remArgs.append("--help")
    if args.verbose:
        remArgs.append("--verbose")

    if args.verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    modArg = args.module
    clName = None
    if ":" in modArg:
        if modArg.count(":") > 1:
            raise RuntimeError("More than one : in module argument '{0}'".format(modArg))
        modArg, clName = tuple(modArg.split(":"))
    ## modArg should now be a python module name or path now
    spec = None
    try:
        spec = importlib.util.find_spec(modArg)
    except:
        pass
    if spec: ## is a module
        if not clName:
            clName = spec.split(".")[-1]
    elif os.path.exists(modArg):
        modName = os.path.splitext(os.path.basename(modArg))[0]
        if not clName:
            clName = modName
        spec = importlib.util.spec_from_file_location(modName, modArg)
    else:
        raise RuntimeError("Module argument '{0}' is neither an importable module or an existing path".format(modArg))
    ## do the import
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, clName):
        raise RuntimeError("Module {0} has no class with name {1}, please specify it using --module=modulenameOrPath:classname")
    modCls = getattr(mod, clName)
    modInst = modCls(remArgs)
    modInst.run()

if __name__ == "__main__":
    main()
