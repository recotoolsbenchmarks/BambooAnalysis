"""
math/builtins-like module for tree operations
"""
## used as a namespace, so avoid filling it with lower-level objects
from . import treeoperations as _to
from . import treeproxies as _tp
import sys
import os.path
pkgRoot = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
import ROOT
ROOT.gROOT.ProcessLine('#include "Math/VectorUtil.h"')
pkgPrefix = os.path.dirname(os.path.dirname(os.path.dirname(pkgRoot)))
ROOT.gInterpreter.AddIncludePath(os.path.join(pkgPrefix, "include", "site", "python{0.major}.{0.minor}".format(sys.version_info), "bamboo"))
ROOT.gSystem.Load(os.path.join(pkgRoot, "libBinnedValues"))
for fname in ("range.h", "jmesystematics.h", "scalefactors.h"):
    ROOT.gROOT.ProcessLine('#include "{}"'.format(fname))
del ROOT,os,sys

## simple type support
def c_bool(arg):
    if arg:
        return _tp.makeConst("true", _tp.boolType)
    else:
        return _tp.makeConst("false", _tp.boolType)
def c_int(num):
    return _tp.makeConst(num, "Int_t")
def c_float(num):
    return _tp.makeConst(num, "Float_t")

## boolean logic
def NOT(sth):
    return _to.MathOp("not", sth, outType=_tp.boolType).result
def AND(*args):
    return _to.MathOp("and", *args, outType=_tp.boolType).result
def OR(*args):
    return _to.MathOp("or", *args, outType=_tp.boolType).result

def switch(test, trueBranch, falseBranch):
    assert trueBranch._typeName == falseBranch._typeName
    return _to.MathOp("switch", test, trueBranch, falseBranch, outType=trueBranch._typeName)
def extMethod(name):
    return _tp.MethodProxy(name) ## TODO somehow take care of includes as well
def extVar(typeName, name):
    return _to.ExtVar(typeName, name).result
def construct(typeName, args):
    return _to.Construct(typeName, args).result
def initList(typeName, elmName, elms):
    return _to.InitList(typeName, elmName, elms).result
def define(typeName, definition):
    return _to.DefinedVar(typeName, definition).result

## math
def abs(sth):
    return _to.MathOp("abs", sth, outType="Float_t").result
def sign(sth):
    return switch(sth!=0., sth/abs(sth), 0.)
def sum(*args, **kwargs):
    return _to.MathOp("add", *args, outType=kwargs.pop("outType", "Float_t")).result
def product(*args):
    return _to.MathOp("multiply", *args, outType="Float_t").result
def log(sth):
    return _to.MathOp("log", sth, outType="Float_t").result
def log10(sth):
    return _to.MathOp("log10", sth, outType="Float_t").result
def max(a1,a2):
    return _to.MathOp("max", a1, a2, outType="Float_t").result
def min(a1,a2):
    return _to.MathOp("min", a1, a2, outType="Float_t").result
def in_range(low, arg, up):
    return _to.AND(arg > low, arg < up)

## Kinematics and helpers
def invariant_mass(*args):
    return extMethod("ROOT::Math::VectorUtil::InvariantMass")(*args)
def invariant_mass_squared(*args):
    return extMethod("ROOT::Math::VectorUtil::InvariantMass2")(*args)
def deltaPhi(a1, a2):
    return extMethod("Kinematics::deltaPhi")(a1, a2)
def deltaR(a1, a2):
    return extMethod("Kinematics::deltaR")(a1, a2)
def signedDeltaPhi(a1, a2):
    return extMethod("Kinematics::signedDeltaPhi")(a1, a2)
def signedDeltaEta(a1, a2):
    return extMethod("Kinematics::signedDeltaEta")(a1, a2)


## range operations
def rng_len(sth):
    return sth.__len__() ## __builtins__.len checks it is an integer

def rng_sum(rng, fun=lambda x : x, start=c_float(0.)):
    return _to.Reduce(rng, start, ( lambda fn : (
        lambda res, elm : res+fn(elm)
        ) )(fun) ).result

def rng_count(rng, pred=lambda x : c_bool(True)): ## specialised version of sum, for convenience
    return _to.Reduce(rng, c_int(0), ( lambda prd : (
        lambda res, elm : res+switch(prd(elm), c_int(1), c_int(0))
        ) )(pred) ).result

def rng_product(rng, fun=lambda x : x):
    return _to.Reduce(rng, c_float(1.), ( lambda fn : (
        lambda res, elm : res*fn(elm)
        ) )(fun) ).result

def rng_max(rng, fun=lambda x : x):
    return _to.Reduce(rng, c_float(float("-inf")), ( lambda fn : (
        lambda res, elm : max(res, fn(elm))
        ) )(fun) ).result
def rng_min(rng, fun=lambda x : x):
    return _to.Reduce(rng, c_float(float("+inf")), ( lambda fn : (
        lambda res, elm : min(res, fn(elm))
        ) )(fun) ).result

def rng_max_element_by(rng, fun=lambda elm : elm):
    return rng._base[_to.Reduce(rng, _tp.makeConst(0, _to.SizeType), ( lambda fn,rn : (
        lambda ibest, elm : switch(fn(elm) > fn(rn._base[ibest]), elm._idx.result, ibest)
        ) )(fun,rng) ).result]
def rng_min_element_by(rng, fun=lambda elm : elm):
    return rng._base[_to.Reduce(rng, _tp.makeConst(0, _to.SizeType), ( lambda fn,rn : (
        lambda ibest, elm : switch(fn(elm) < fn(rn._base[ibest]), elm._idx.result, ibest)
        ) )(fun,rng) ).result]

## early-exit algorithms
def rng_any(rng, pred=lambda elm : elm):
    return rng.__len__() != _to.Next(rng, pred)
def rng_find(rng, pred=lambda elm : _tp.makeConst("true", boolType)):
    return _to.Next(rng, pred).result

## range-to-range: selection, combinatorics, systematic variations
def select(rng, pred=lambda elm : True):
    return _tp.SelectionProxy(_to.Select(rng, pred), rng)


def addKinematicVariation(vrng, key, modif=( lambda elm : elm.p4 ), pred=( lambda mp4,elm : True )):
    modif = _to.KinematicVariation(vrng.nominal, modif, pred).result
    modif.itemType = vrng._varItemType
    vrng.register(key, modif)

def combine(rng, N=None, pred=lambda combo : True):
    ## TODO figure out if 'rng' is one range, or N
    ## - if one -> make combinations of similar particles (e.g. Z->l+l-, W->jj)
    ## - if N   -> make combinations of different particles (e.g. t->lj)
    return NotImplemented ## store tuples of (smart) refs if they pass the cut