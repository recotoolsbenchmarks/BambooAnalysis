## used as a namespace, so avoid filling it with lower-level objects
from . import treeoperations as _to
from . import treeproxies as _tp

def _load_extensions():
    """Add extension libraries and necessary header files to the ROOT interpreter"""
    import sys
    import os.path
    pkgRoot = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    from cppyy import gbl
    gbl.gROOT.ProcessLine('#include "Math/VectorUtil.h"')
    instPrefix = os.path.dirname(os.path.dirname(os.path.dirname(pkgRoot)))
    instInclude = os.path.join(instPrefix, "include", "site", "python{0.major}.{0.minor}".format(sys.version_info), "bamboo")
    if os.path.isdir(instInclude): ## installed mode
        gbl.gInterpreter.AddIncludePath(instInclude)
        libDir = pkgRoot
    else: ## non-installed mode
        libDir = os.path.join(pkgRoot, "build", "lib")
        if not os.path.isdir(libDir):
            raise RuntimeError("No directory {0} so running in local mode, but then build/lib need to be present. Did you run 'python setup.py build'?")
        gbl.gInterpreter.AddIncludePath(os.path.join(pkgRoot, "build", "include"))
        gbl.gInterpreter.AddIncludePath(os.path.join(pkgRoot, "cpp"))
    gbl.gSystem.Load(os.path.join(libDir, "libBinnedValues"))
    gbl.gSystem.Load(os.path.join(libDir, "libBambooLumiMask"))
    for fname in ("range.h", "jmesystematics.h", "scalefactors.h", "LumiMask.h"):
        gbl.gROOT.ProcessLine('#include "{}"'.format(fname))
_load_extensions()

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
    if len(args) == 1:
        return _tp.makeProxy(_tp.boolType, _to.adaptArg(args[0]))
    else:
        return _to.MathOp("and", *args, outType=_tp.boolType).result
def OR(*args):
    if len(args) == 1:
        return _tp.makeProxy(_tp.boolType, _to.adaptArg(args[0]))
    else:
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
    """ Add a container with modified kinematics (for jet systematics)

    :param vrng: base object range (kinematic variation-enabled, see :py:class:`bamboo.treeproxies.Variations`)
    :param key: registered name of the variation
    :param modif: transformation to apply, object->modified-p4
    :param pred: selection to apply after performing the transformation
    """
    modif = _to.KinematicVariation(vrng.nominal, modif, pred).result
    modif.itemType = vrng._varItemType
    vrng.register(key, modif)

def combine(rng, N=None, pred=lambda *parts : True, sameIdxPred=lambda i1,i2: i1 < i2):
    """ Create N-particle combination from one or several ranges

    :param rng: range (or iterable of ranges) with basic objects to combine
    :param N: number of objects to combine (at least 2), in case of multiple ranges it does not need to be given
    (``len(rng)`` will be taken; if specified they should match)
    :param pred: selection to apply to candidates
    :param sameIdxPred: additional selection for objects from the same base container (by index).
        The default avoids duplicates by keeping the indices sorted
        (so the first object in the container will also be the first of the combination),
        but in some cases one may want to simply test inequality.

    :Example:

    >>> osdimu = op.combine(t.Muon, N=2, pred=lambda mu1,mu2 : mu1.charge != mu2.charge)
    >>> oselmu = op.combine((t.Electron, t.Muon), pred=lambda el,mu : el.charge != mu.charge)

    .. note::

        currently only N=2 is supported
    """
    if not hasattr(rng, "__iter__"):
        rng = (rng,)
    elif N is None:
        N = len(rng)
    elif N <= 1:
        raise RuntimeError("Can only make combinations of more than one")
    if len(rng) != N and len(rng) != 1:
        raise RuntimeError("If N(={0:d}) input ranges are given, only N-combinations can be made, not {1:d}".format(len(rng), N))
    return _to.Combine(N, rng, pred, sameIdxPred=sameIdxPred).result ## only implemented for 2 (same or different container)
