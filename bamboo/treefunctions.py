## used as a namespace, so avoid filling it with lower-level objects
from . import treeoperations as _to
from . import treeproxies as _tp

def _load_extensions():
    """Add extension libraries and necessary header files to the ROOT interpreter"""
    import sys
    import pkg_resources
    import os.path
    pkgRoot = pkg_resources.get_distribution("bamboo").location
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
            raise RuntimeError("No directory {0} so running in local mode, but then build/lib need to be present. Did you run 'python setup.py build'?".format(libDir))
        gbl.gInterpreter.AddIncludePath(os.path.join(pkgRoot, "build", "include"))
        gbl.gInterpreter.AddIncludePath(os.path.join(pkgRoot, "cpp"))
    gbl.gSystem.Load(os.path.join(libDir, "libBinnedValues"))
    gbl.gSystem.Load(os.path.join(libDir, "libBambooLumiMask"))
    gbl.gSystem.Load(os.path.join(libDir, "libJMEObjects"))
    gbl.gSystem.Load(os.path.join(libDir, "libRoccoR"))
    ## TODO combine into libBamboo, and "bamboo.h"?
    for fname in ("bamboohelpers.h", "range.h", "scalefactors.h", "LumiMask.h", "JMESystematicsCalculator.h", "RochesterCorrectionCalculator.h"):
        gbl.gROOT.ProcessLine('#include "{}"'.format(fname))
    getattr(gbl, "JMESystematicsCalculator::result_t") ## trigger dictionary generation
_load_extensions()

## simple type support
def c_bool(arg):
    """ Construct a boolean constant """
    if arg:
        return _tp.makeProxy(_tp.boolType, _to.Const(_tp.boolType, "true"))
    else:
        return _tp.makeProxy(_tp.boolType, _to.Const(_tp.boolType, "false"))
def c_int(num):
    """ Construct an integer number constant """
    return _tp.makeConst(num, "Int_t")
def c_float(num):
    """ Construct a floating-point number constant """
    return _tp.makeConst(num, "Float_t")

## boolean logic
def NOT(sth):
    """ Logical NOT """
    return _to.MathOp("not", sth, outType=_tp.boolType).result
def AND(*args):
    """ Logical AND """
    if len(args) == 1:
        return _tp.makeProxy(_tp.boolType, _to.adaptArg(args[0]))
    else:
        return _to.MathOp("and", *args, outType=_tp.boolType).result
def OR(*args):
    """ Logical OR """
    if len(args) == 1:
        return _tp.makeProxy(_tp.boolType, _to.adaptArg(args[0]))
    else:
        return _to.MathOp("or", *args, outType=_tp.boolType).result

def switch(test, trueBranch, falseBranch):
    """ Pick one or another value, based on a third one (ternary operator in C++)

    :Example:

    >>> op.switch(runOnMC, mySF, 1.) ## incomplete pseudocode
    """
    # assert trueBranch._typeName == falseBranch._typeName
    return _to.MathOp("switch", test, trueBranch, falseBranch, outType=trueBranch._typeName)
def extMethod(name):
    return _tp.MethodProxy(name) ## TODO somehow take care of includes as well
def extVar(typeName, name):
    return _to.ExtVar(typeName, name).result
def construct(typeName, args):
    if not hasattr(args, "__iter__"):
        args = (args,)
    return _to.Construct(typeName, args).result
def static_cast(typeName, arg):
    return _tp.makeProxy(typeName, _to.CallMethod("static_cast<{0}>".format(typeName), (arg,), getFromRoot=False))
def initList(typeName, elmName, elms):
    return _to.InitList(typeName, elms, elmType=elmName).result
def define(typeName, definition, nameHint=None):
    return _to.DefinedVar(typeName, definition, nameHint=nameHint).result

## math
def abs(sth):
    """ Return the absolute value

    :Example:

    >>> op.abs(t.Muon[0].p4.Eta())
    """
    return _to.MathOp("abs", sth, outType=sth._typeName).result
def sign(sth):
    """ Return the sign of a number

    :Example:

    >>> op.sign(t.Muon[0].p4.Eta())
    """
    return switch(sth!=0., sth/abs(sth), 0.)
def sum(*args, **kwargs):
    """ Return the sum of the arguments

    :Example:

    >>> op.sum(t.Muon[0].p4.Eta(), t.Muon[1].p4.Eta())
    """
    return _to.MathOp("add", *args, outType=kwargs.pop("outType",
        _tp.floatType if _tp._hasFloat(*args) else _tp.intType)).result
def product(*args):
    """ Return the product of the arguments

    :Example:

    >>> op.product(t.Muon[0].p4.Eta(), t.Muon[1].p4.Eta())
    """
    return _to.MathOp("multiply", *args,
            outType=(_tp.floatType if _tp._hasFloat(*args) else _tp.intType)).result
def sqrt(sth):
    """ Return the square root of a number

    :Example:

    >>> m1, m2 = t.Muon[0].p4, t.Muon[1].p4
    >>> m12dR = op.sqrt( op.pow(m1.Eta()-m2.Eta(), 2) + op.pow(m1.Phi()-m2.Phi(), 2) )
    """
    return _to.MathOp("sqrt", sth, outType=sth._typeName).result
def pow(base, exp):
    """ Return a power of a number

    :Example:

    >>> m1, m2 = t.Muon[0].p4, t.Muon[1].p4
    >>> m12dR = op.sqrt( op.pow(m1.Eta()-m2.Eta(), 2) + op.pow(m1.Phi()-m2.Phi(), 2) )
    """
    return _to.MathOp("pow", base, exp, outType=base._typeName).result
def exp(sth):
    """ Return the exponential of a number

    :Example:

    >>> op.exp(op.abs(t.Muon[0].p4.Eta()))
    """
    return _to.MathOp("exp", sth).result
def log(sth):
    """ Return the natural logarithm of a number

    :Example:

    >>> op.log(t.Muon[0].p4.Pt())
    """
    return _to.MathOp("log", sth).result
def log10(sth):
    """ Return the base-10 logarithm of a number

    :Example:

    >>> op.log10(t.Muon[0].p4.Pt())
    """
    return _to.MathOp("log10", sth).result

def sin(sth):
    """ Return the sine of a number

    :Example:

    >>> op.sin(t.Muon[0].p4.Phi())
    """
    return _to.MathOp("sin", sth).result

def cos(sth):
    """ Return the cosine of a number

    :Example:

    >>> op.cos(t.Muon[0].p4.Phi())
    """
    return _to.MathOp("cos", sth).result

def tan(sth):
    """ Return the tangent of a number

    :Example:

    >>> op.tan(t.Muon[0].p4.Phi())
    """
    return _to.MathOp("tan", sth).result

def asin(sth):
    """ Return the arcsine of a number

    :Example:

    >>> op.asin(op.c_float(3.1415))
    """
    return _to.MathOp("asin", sth).result

def acos(sth):
    """ Return the arccosine of a number

    :Example:

    >>> op.ascos(op.c_float(3.1415))
    """
    return _to.MathOp("acos", sth).result

def atan(sth):
    """ Return the arctangent of a number

    :Example:

    >>> op.atan(op.c_float(3.1415))
    """
    return _to.MathOp("atan", sth).result

def max(a1,a2):
    """ Return the maximum of two numbers

    :Example:

    >>> op.max(op.abs(t.Muon[0].p4.Eta()), op.abs(t.Muon[1].p4.Eta()))
    """
    return _to.MathOp("max", a1, a2,
            outType=(_tp.floatType if _tp._hasFloat(a1, a2) else _tp.intType)).result
def min(a1,a2):
    """ Return the minimum of two numbers

    :Example:

    >>> op.min(op.abs(t.Muon[0].p4.Eta()), op.abs(t.Muon[1].p4.Eta()))
    """
    return _to.MathOp("min", a1, a2,
            outType=(_tp.floatType if _tp._hasFloat(a1, a2) else _tp.intType)).result
def in_range(low, arg, up):
    """ Check if a value is inside a range (boundaries excluded)

    :Example:

    >>> op.in_range(10., t.Muon[0].p4.Pt(), 20.)
    """
    return extMethod("rdfhelpers::in_range")(*(_to.adaptArg(iarg, typeHint=_tp.floatType) for iarg in (low, arg, up)))

## Kinematics and helpers
def withMass(arg, massVal):
    """ Construct a Lorentz vector with given mass (taking the other components from the input)

    :Example:

    >>> pW = withMass((j1.p4+j2.p4), 80.4)
    """
    return extMethod("rdfhelpers::withMass")(arg, _to.adaptArg(massVal, typeHint=_tp.floatType))
def invariant_mass(*args):
    """ Calculate the invariant mass of the arguments

    :Example:

    >>> mElEl = op.invariant_mass(t.Electron[0].p4, t.Electron[1].p4)

    .. note::

        Unlike in the example above, :py:meth:`bamboo.treefunctions.combine` should be used to make N-particle combinations in most practical cases
    """
    if len(args) == 0:
        raise RuntimeError("Need at least one argument to calculate invariant mass")
    elif len(args) == 1:
        return args[0].M()
    elif len(args) == 2:
        return extMethod("ROOT::Math::VectorUtil::InvariantMass")(*args)
    else:
        return sum(*args, outType=args[0]._typeName).M()
def invariant_mass_squared(*args):
    """ Calculate the squared invariant mass of the arguments using ``ROOT::Math::VectorUtil::InvariantMass2``

    :Example:

    >>> m2ElEl = op.invariant_mass2(t.Electron[0].p4, t.Electron[1].p4)
    """
    if len(args) == 0:
        raise RuntimeError("Need at least one argument to calculate invariant mass")
    elif len(args) == 1:
        return args[0].M2()
    elif len(args) == 2:
        return extMethod("ROOT::Math::VectorUtil::InvariantMass2")(*args)
    else:
        return sum(*args, outType=args[0]._typeName).M2()
def deltaPhi(a1, a2):
    """ Calculate the difference in azimutal angles (using ``ROOT::Math::VectorUtil::DeltaPhi``)

    :Example:

    >>> elelDphi = op.deltaPhi(t.Electron[0].p4, t.Electron[1].p4)
    """
    return extMethod("ROOT::Math::VectorUtil::DeltaPhi")(a1, a2)
def deltaR(a1, a2):
    """ Calculate the Delta R distance (using ``ROOT::Math::VectorUtil::DeltaR``)

    :Example:

    >>> elelDR = op.deltaR(t.Electron[0].p4, t.Electron[1].p4)
    """
    return extMethod("ROOT::Math::VectorUtil::DeltaR")(a1, a2)

## range operations
def rng_len(sth):
    """ Get the number of elements in a range

    :param rng: input range

    :Example:

    >>> nElectrons = op.rng_len(t.Electron)
    """
    return sth.__len__() ## __builtins__.len checks it is an integer

def rng_sum(rng, fun=lambda x : x, start=c_float(0.)):
    """ Sum the values of a function over a range

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> totalMuCharge = op.rng_sum(t.Muon, lambda mu : mu.charge)
    """
    return _to.Reduce.fromRngFun(rng, start, ( lambda fn : (
        lambda res, elm : res+fn(elm)
        ) )(fun) )

def rng_count(rng, pred=lambda x : c_bool(True)): ## specialised version of sum, for convenience
    """ Count the number of elements passing a selection

    :param rng: input range
    :param pred: selection predicate (a callable that takes an element of the range and returns a boolean)

    :Example:

    >>> nCentralMu = op.rng_count(t.Muon, lambda mu : op.abs(mu.p4.Eta() < 2.4))
    """
    return _to.Reduce.fromRngFun(rng, c_int(0), ( lambda prd : (
        lambda res, elm : res+switch(prd(elm), c_int(1), c_int(0))
        ) )(pred) )

def rng_product(rng, fun=lambda x : x):
    """ Calculate the production of a function over a range

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> overallMuChargeSign = op.rng_product(t.Muon, lambda mu : mu.charge)
    """
    return _to.Reduce.fromRngFun(rng, c_float(1.), ( lambda fn : (
        lambda res, elm : res*fn(elm)
        ) )(fun) )

def rng_max(rng, fun=lambda x : x):
    """ Find the highest value of a function in a range

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> mostForwardMuEta = op.rng_max(t.Muon. lambda mu : op.abs(mu.p4.Eta()))
    """
    return _to.Reduce.fromRngFun(rng, c_float(float("-inf")), ( lambda fn : (
        lambda res, elm : extMethod("std::max")(res, fn(elm))
        ) )(fun) )
def rng_min(rng, fun=lambda x : x):
    """ Find the lowest value of a function in a range

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> mostCentralMuEta = op.rng_min(t.Muon. lambda mu : op.abs(mu.p4.Eta()))
    """
    return _to.Reduce.fromRngFun(rng, c_float(float("+inf")), ( lambda fn : (
        lambda res, elm : extMethod("std::min")(res, fn(elm))
        ) )(fun) )

def rng_max_element_by(rng, fun=lambda elm : elm):
    """ Find the element for which the value of a function is maximal

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> mostForwardMu = op.rng_max_element_by(t.Muon. lambda mu : op.abs(mu.p4.Eta()))
    """
    return rng._base[_to.Reduce.fromRngFun(rng,
        construct("std::pair<{0},{1}>".format(_tp.SizeType,_tp.floatType), (c_int(-1), c_float(float("-inf")))),
        ( lambda fn : ( lambda ibest, elm : extMethod("rdfhelpers::maxPairBySecond")(ibest, elm._idx.result, fn(elm)) ) )(fun)).first]

def rng_min_element_by(rng, fun=lambda elm : elm):
    """ Find the element for which the value of a function is minimal

    :param rng: input range
    :param fun: function whose value should be used (a callable that takes an element of the range and returns a number)

    :Example:

    >>> mostCentralMu = op.rng_min_element_by(t.Muon. lambda mu : op.abs(mu.p4.Eta()))
    """
    return rng._base[_to.Reduce.fromRngFun(rng,
        construct("std::pair<{0},{1}>".format(_tp.SizeType,_tp.floatType), (c_int(-1), c_float(float("+inf")))),
        ( lambda fn : ( lambda ibest, elm : extMethod("rdfhelpers::minPairBySecond")(ibest, elm._idx.result, fn(elm)) ) )(fun)).first]

## early-exit algorithms
def rng_any(rng, pred=lambda elm : elm):
    """ Test if any item in a range passes a selection

    :param rng: input range
    :param pred: selection predicate (a callable that takes an element of the range and returns a boolean)

    :Example:

    >>> hasCentralMu = op.rng_any(t.Muon. lambda mu : op.abs(mu.p4.Eta()) < 2.4)
    """
    return  _tp.makeConst(-1, _to.SizeType) != _to.Next.fromRngFun(rng, pred)._idx
def rng_find(rng, pred=lambda elm : _tp.makeConst("true", boolType)):
    """ Find the first item in a range that passes a selection

    :param rng: input range
    :param pred: selection predicate (a callable that takes an element of the range and returns a boolean)

    :Example:

    >>> leadCentralMu = op.rng_find(t.Muon, lambda mu : op.abs(mu.p4.Eta()) < 2.4)
    """
    return _to.Next.fromRngFun(rng, pred)

## range-to-range: selection, combinatorics, systematic variations
def select(rng, pred=lambda elm : True):
    """ Select elements from the range that pass a cut

    :param rng: input range
    :param pred: selection predicate (a callable that takes an element of the range and returns a boolean)

    :Example:

    >>> centralMuons = op.select(t.Muon, lambda mu : op.abs(mu.p4.Eta()) < 2.4)
    """
    return _to.Select.fromRngFun(rng, pred)

def sort(rng, fun=lambda elm : 0.):
    """ Sort the range (ascendingly) by the value of a function applied on each element

    :param rng: input range
    :param fun: function by whose value the elements should be sorted

    :Example:

    >>> muonsByCentrality = op.sort(t.Muon, lambda mu : op.abs(mu.p4.Eta()))
    """
    return _to.Sort.fromRngFun(rng, fun)

def map(rng, fun, valueType=None):
    """ Create a list of derived values for a collection

    This is useful for storing a derived quantity each item of a collection on a skim,
    and also for filling a histogram for each entry in a collection.

    :param rng: input range
    :param fun: function to calculate derived values
    :param valueType: stored return type (optional, ``fun(rng[i])`` should be convertible to this type)

    :Example:

    >>> muon_absEta = op.map(t.Muon, lambda mu : op.abs(mu.p4.Eta()))
    """
    return _to.Map.fromRngFun(rng, fun, typeName=valueType)

def rng_pickRandom(rng, seed=0):
    """ Pick a random element from a range

    :param rng: range to pick an element from
    :param seed: seed for the random generator

    .. caution::

        empty placeholder, to be implemented
    """
    return rng[_to.PseudoRandom(0, rng_len(rng), seed, isIntegral=True)]

def combine(rng, N=None, pred=lambda *parts : True, sameIdxPred=lambda i1,i2: i1 < i2):
    """ Create N-particle combination from one or several ranges

    :param rng: range (or iterable of ranges) with basic objects to combine
    :param N: number of objects to combine (at least 2), in case of multiple ranges it does not need to be given (``len(rng)`` will be taken; if specified they should match)
    :param pred: selection to apply to candidates (a callable that takes the constituents and returns a boolean)
    :param sameIdxPred: additional selection for objects from the same base container (by index; a callable that takes the indices and returns a boolean). The default avoids duplicates by keeping the indices sorted (so the first object in the container will also be the first of the combination), but in some cases one may want to simply test inequality.

    :Example:

    >>> osdimu = op.combine(t.Muon, N=2, pred=lambda mu1,mu2 : mu1.charge != mu2.charge)
    >>> firstosdimu = osdimu[0]
    >>> firstosdimu_Mll = op.invariant_mass(firstosdimu[0].p4, firstosdimu[1].p4)
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
    return _to.Combine.fromRngFun(N, rng, pred, sameIdxPred=sameIdxPred) ## only implemented for 2 (same or different container)
