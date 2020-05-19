"""
Object representation of operations on TTree branches

The aim is to provide provide sufficiently flexible and complete foundation
for the development of efficient histogram-filling programs
through the use of python wrappers (see e.g. treeproxies).
"""

from itertools import chain, repeat, combinations, count, tee
from contextlib import contextmanager
import logging
logger = logging.getLogger(__name__)

class TupleOpCache:
    __slots__ = ("hash", "repr")
    def __init__(self):
        self.hash = None
        self.repr = None
    def __bool__(self):
        return self.hash is not None or self.repr is not None

class TupleOp:
    """ Interface & base class for operations on leafs and resulting objects / values

    Instances should be defined once, and assumed immutable by all observers
    (they should only ever be modified just after construction, preferably by the owner).
    Since a value-based hash (and repr) is cached on first use, violating this rule
    might lead to serious bugs. In case of doubt the clone method can be used to
    obtain an independent copy.
    Subclasses should define a result property and clone, _repr, _eq, and optionally deps methods.
    """
    __slots__ = ("_cache", "canDefine") ## this means all deriving classes need to define __slots__ (for performance)
    def __init__(self):
        self._cache = TupleOpCache()
        self.canDefine = True
    def clone(self, memo=None, select=lambda nd : True):
        """ Create an independent copy (with empty repr/hash cache) of the (sub)expression """
        if memo is None: ## top-level, construct the dictionary
            memo = dict()
        if id(self) in memo:
            return memo[id(self)]
        else:
            cp = self._clone(memo, select)
            cln = cp if cp is not None else self
            memo[id(self)] = cln
            return cln
    def _clone(self, memo, select): ## simple version, call clone of attributes without worrying about memo
        """ Implementation of clone - to be overridden by all subclasses (memo is dealt with by clone, so simply construct, calling .clone(memo=memo, select=select) on TupleOp attributes. Return None if no clone is needed (based on select) """
        if select(self):
            return self.__class__()
    def deps(self, defCache=None, select=(lambda x : True), includeLocal=False):
        """ Dependent TupleOps iterator """
        yield from []
    @property
    def result(self):
        """ Proxy to the result of this (sub)expression """
        pass
    ## subclasses should define at least _clone, _repr, and _eq (value-based)
    def __repr__(self):
        """ String representation (used for hash, and lazily cached) """
        if self._cache.repr is None:
            self._cache.repr = self._repr()
        return self._cache.repr
    def _repr(self):
        """ __repr__ implementation - to be overridden by all subclasses (caching is in top-level __repr__) """
        return "TupleOp()"
    def __hash__(self):
        """ Value-based hash (lazily cached) """
        if self._cache.hash is None:
            self._cache.hash = self._hash()
        return self._cache.hash
    def _hash(self):
        """ __hash__ implementation - subclasses please override  (caching is in top-level __hash__) """
        return hash(self.__repr__())
    def __eq__(self, other):
        """ Identity or value-based equality comparison (same object and unequal should be fast) """
        # _eq may end up being quite expensive, but should almost never be called
        return id(self) == id(other) or ( self.__hash__() == hash(other) and self.__class__ == other.__class__ and self._eq(other) )
    def _eq(self, other):
        """ value-based __eq__ implementation - to be overridden by all subclasses (protects against hash collisions; hash and class are checked to be equal already) """
        return True
    def get_cppStr(self, defCache=None):
        """ Interface method: generate a C++ string from the operation

        :param defCache: cache with defined operations and symbols; a minimal implementation is :py:class:`~bamboo.treeoperations.CppStrRedir` (which does not do any caching)
        """
        pass

## implementations are split out, see treeproxies
class TupleBaseProxy:
    """
    Interface & base class for proxies
    """
    def __init__(self, typeName, parent=None):
        self._typeName = typeName
        self._parent = parent
    @property
    def op(self):
        if self._parent is None:
            raise ValueError("Cannot get operation for {0!r}, abstract base class / empty parent".format(self))
        return self._parent

class CppStrRedir:
    """ Expression cache interface. Default implementation: no caching """
    def __init__(self):
        self._iFun = 0
    def __call__(self, arg):
        return arg.get_cppStr(defCache=self)
    def symbol(self, decl):
        """
        Define (or get) a new C++ symbol for the declaration

        decl should contain the code, with <<name>> where the name should go.  Returns the unique name
        """
        print("WARNING: should add defined symbol for '{0}' but that's not supported".format(decl))
    def _getColName(self, op):
        return None

cppNoRedir = CppStrRedir()

class ForwardingOp(TupleOp):
    """ Transparent wrapper (base for marking parts of the tree, e.g. things with systematic variations) """
    __slots__ = ("wrapped",)
    def __init__(self, wrapped, canDefine=None):
        super(ForwardingOp, self).__init__()
        self.wrapped = wrapped
        self.canDefine = canDefine if canDefine is not None else self.wrapped.canDefine
    def _clone(self, memo, select):
        clWr = self.wrapped.clone(memo=memo, select=select)
        if select(self) or id(clWr) != id(self.wrapped):
            return self.__class__(clWr, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if select(self.wrapped):
            yield self.wrapped
        yield from self.wrapped.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        wrapRes = self.wrapped.result
        wrapRes._parent = self
        return wrapRes
    def _repr(self):
        return "{0}({1!r})".format(self.__class__.__name__, self.wrapped)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.wrapped)))
    def _eq(self, other):
        return self.wrapped == other.wrapped
    def get_cppStr(self, defCache=cppNoRedir):
        return defCache(self.wrapped)

SizeType = "std::size_t"

class Const(TupleOp):
    """ Hard-coded number (or expression) """
    __slots__ = ("typeName", "value")
    def __init__(self, typeName, value):
        super(Const, self).__init__()
        self.typeName = typeName
        self.value = value
    def _clone(self, memo, select):
        if select(self):
            return self.__class__(self.typeName, self.value)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.value)
    def _hash(self):
        return hash(self._repr())
    def _eq(self, other):
        return self.typeName == other.typeName and self.value == other.value
    # backends
    def get_cppStr(self, defCache=None):
        try:
            if abs(self.value) == float("inf"):
                return "std::numeric_limits<{0}>::{mnmx}()".format(self.typeName, mnmx=("min" if self.value < 0. else "max"))
        except:
            pass
        return str(self.value) ## should maybe be type-aware...

class GetColumn(TupleOp):
    """ Get a column value """
    __slots__ = ("typeName", "name")
    def __init__(self, typeName, name):
        super(GetColumn, self).__init__()
        self.typeName = typeName
        self.name = name
    def _clone(self, memo, select):
        if select(self):
            return self.__class__(self.typeName, self.name)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.name)
    def _hash(self):
        return hash(self._repr())
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name
    def get_cppStr(self, defCache=None):
        return self.name

class GetArrayColumn(TupleOp):
    """ Get the number from a leaf """
    __slots__ = ("typeName", "name", "length")
    def __init__(self, typeName, name, length):
        super(GetArrayColumn, self).__init__()
        self.typeName = typeName
        self.name = name
        self.length = length
    def _clone(self, memo, select):
        lnCl = self.length.clone(memo=memo, select=select)
        if select(self) or id(lnCl) != id(self.length):
            return self.__class__(self.typeName, self.name, lnCl)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            if select(self.length):
                yield self.length
            yield from self.length.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self, length=makeProxy(SizeType, self.length))
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.typeName, self.name, self.length)
    def _hash(self):
        return hash(self._repr())
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name and self.length == other.length
    def get_cppStr(self, defCache=None):
        return self.name

## helper
def adaptArg(arg, typeHint=None):
    if isinstance(arg, TupleBaseProxy):
        return arg.op
    elif isinstance(arg, TupleOp):
        return arg
    elif typeHint is not None:
        if str(arg) == arg: ## string, needs quote
            return Const(typeHint, '"{}"'.format(arg))
        elif isinstance(arg, bool):
            return Const(typeHint, "true" if arg else "false")
        else:
            return Const(typeHint, arg)
    else:
        raise ValueError("Should get some kind of type hint")

mathOpFuns_cppStr = {
      "add"      : lambda cppStr,*args : "( {0} )".format(" + ".join(cppStr(arg) for arg in args))
    , "multiply" : lambda cppStr,*args : "( {0} )".format(" * ".join(cppStr(arg) for arg in args))
    , "subtract" : lambda cppStr,a1,a2 : "( {0} - {1} )".format(cppStr(a1), cppStr(a2))
    , "divide"   : lambda cppStr,a1,a2 : "( {0} / {1} )".format(cppStr(a1), cppStr(a2))
    , "floatdiv"  : lambda cppStr,a1,a2 : "( 1.*{0} / {1} )".format(cppStr(a1), cppStr(a2))
    , "mod" : lambda cppStr,a1,a2 : "( {0} % {1} )".format(cppStr(a1), cppStr(a2))
    , "neg"  : lambda cppStr,a : "( -{0} )".format(cppStr(a))
    #
    , "lt" : lambda cppStr,a1,a2 : "( {0} <  {1} )".format(cppStr(a1), cppStr(a2))
    , "le" : lambda cppStr,a1,a2 : "( {0} <= {1} )".format(cppStr(a1), cppStr(a2))
    , "eq" : lambda cppStr,a1,a2 : "( {0} == {1} )".format(cppStr(a1), cppStr(a2))
    , "ne" : lambda cppStr,a1,a2 : "( {0} != {1} )".format(cppStr(a1), cppStr(a2))
    , "gt" : lambda cppStr,a1,a2 : "( {0} >  {1} )".format(cppStr(a1), cppStr(a2))
    , "ge" : lambda cppStr,a1,a2 : "( {0} >= {1} )".format(cppStr(a1), cppStr(a2))
    , "and" : lambda cppStr,*args : "( {0} )".format(" && ".join(cppStr(a) for a in args))
    , "or"  : lambda cppStr,*args : "( {0} )".format(" || ".join(cppStr(a) for a in args))
    , "not" : lambda cppStr,a : "( ! {0} )".format(cppStr(a))
    , "band" : lambda cppStr,*args : "( {0} )".format(" & ".join(cppStr(a) for a in args))
    , "bor"  : lambda cppStr,*args : "( {0} )".format(" | ".join(cppStr(a) for a in args))
    , "bxor"  : lambda cppStr,*args : "( {0} )".format(" ^ ".join(cppStr(a) for a in args))
    , "bnot" : lambda cppStr,a : "( ~ {0} )".format(cppStr(a))
    , "lshift" : lambda cppStr,a1,a2 : "( {0}<<{1} )".format(cppStr(a1), cppStr(a2))
    , "rshift" : lambda cppStr,a1,a2 : "( {0}>>{1} )".format(cppStr(a1), cppStr(a2))
    #
    , "abs"   : lambda cppStr,arg : "std::abs( {0} )".format(cppStr(arg))
    , "sqrt"  : lambda cppStr,arg : "std::sqrt( {0} )".format(cppStr(arg))
    , "pow"   : lambda cppStr,a1,a2 : "std::pow( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "exp"   : lambda cppStr,arg : "std::exp( {0} )".format(cppStr(arg))
    , "log"   : lambda cppStr,arg : "std::log( {0} )".format(cppStr(arg))
    , "log10" : lambda cppStr,arg : "std::log10( {0} )".format(cppStr(arg))
    , "sin"  : lambda cppStr,arg : "std::sin( {0} )".format(cppStr(arg))
    , "cos"  : lambda cppStr,arg : "std::cos( {0} )".format(cppStr(arg))
    , "tan"  : lambda cppStr,arg : "std::tan( {0} )".format(cppStr(arg))
    , "asin"  : lambda cppStr,arg : "std::asin( {0} )".format(cppStr(arg))
    , "acos"  : lambda cppStr,arg : "std::acos( {0} )".format(cppStr(arg))
    , "atan"  : lambda cppStr,arg : "std::atan( {0} )".format(cppStr(arg))
    , "max"   : lambda cppStr,a1,a2 : "std::max( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "min"   : lambda cppStr,a1,a2 : "std::min( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    #
    , "switch" : lambda cppStr,test,trueBr,falseBr : "( ( {0} ) ? ( {1} ) : ( {2} ) )".format(cppStr(test), cppStr(trueBr), cppStr(falseBr))
    }

class MathOp(TupleOp):
    """ Mathematical function N->1, e.g. sin, abs, ( lambda x, y : x*y ) """
    __slots__ = ("outType", "op", "args")
    def __init__(self, op, *args, **kwargs):
        super(MathOp, self).__init__()
        self.outType = kwargs.pop("outType", "Double_t")
        self.op = op
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
        self.canDefine = kwargs.pop("canDefine", all(a.canDefine for a in self.args))
        assert len(kwargs) == 0
    def _clone(self, memo, select):
        argsCl = tuple(a.clone(memo=memo, select=select) for a in self.args)
        if select(self) or any(id(aCl) != id(aOrig) for aCl, aOrig in zip(argsCl, self.args)):
            return self.__class__(self.op, *argsCl, outType=self.outType, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.outType, self)
    def _repr(self):
        return "{0}({1}, {2}, outType={3!r})".format(self.__class__.__name__, self.op, ", ".join(repr(arg) for arg in self.args), self.outType)
    def _hash(self):
        return hash(tuple(chain([ self.__class__.__name__, self.op, self.outType], (hash(a) for a in self.args))))
    def _eq(self, other):
        return self.outType == other.outType and self.op == other.op and self.args == other.args
    def get_cppStr(self, defCache=cppNoRedir):
        return mathOpFuns_cppStr[self.op](defCache, *self.args)

class GetItem(TupleOp):
    """ Get item from array (from function call or from array leaf) """
    __slots__ = ("arg", "typeName", "_index")
    def __init__(self, arg, valueType, index, indexType=SizeType, canDefine=None):
        super(GetItem, self).__init__()
        self.arg = adaptArg(arg)
        self.typeName = valueType
        self._index = adaptArg(index, typeHint=indexType)
        self.canDefine = canDefine if canDefine is not None else ( self.arg.canDefine and self._index.canDefine )
    def _clone(self, memo, select):
        argCl = self.arg.clone(memo=memo, select=select)
        idxCl = self._index.clone(memo=memo, select=select)
        if select(self) or id(argCl) != id(self.arg) or id(self._index) != id(idxCl):
            return self.__class__(argCl, self.typeName, idxCl, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.arg, self._index):
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def index(self):
        from .treeproxies import makeProxy
        return makeProxy(SizeType, self._index)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.arg, self.typeName, self._index)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.arg), self.typeName, hash(self._index)))
    def _eq(self, other):
        return self.arg == other.arg and self.typeName == other.typeName and self._index == other._index
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}[{1}]".format(defCache(self.arg), defCache(self._index))

class Construct(TupleOp):
    __slots__ = ("typeName", "args")
    def __init__(self, typeName, args, canDefine=None):
        super(Construct, self).__init__()
        self.typeName = typeName
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
        self.canDefine = canDefine if canDefine is not None else all(a.canDefine for a in self.args)
    def _clone(self, memo, select):
        argsCl = tuple(a.clone(memo=memo, select=select) for a in self.args)
        if select(self) or any(id(argCl) != id(arg) for argCl, arg in zip(argsCl, self.args)):
            return self.__class__(self.typeName, argsCl, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2})".format(self.__class__.__name__, self.typeName, ", ".join(repr(a) for a in self.args))
    def _hash(self):
        return hash(tuple(chain([ self.__class__.__name__, self.typeName ], (hash(a) for a in self.args))))
    def _eq(self, other):
        return self.typeName == other.typeName and self.args == other.args
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}{{{1}}}".format(self.typeName, ", ".join(defCache(a) for a in self.args))

def guessReturnType(mp):
    if hasattr(mp, "func_doc") and hasattr(mp, "func_name"):
        toks = list(mp.func_doc.split())
        ## left and right strip const * and &
        while toks[-1].rstrip("&") in ("", "const", "static"):
            toks = toks[:-1]
        while toks[0].rstrip("&") in ("", "const", "static"):
            toks = toks[1:]
        while any(tok.endswith("unsigned") for tok in toks):
            iU = next(i for i,tok in enumerate(toks) if tok.endswith("unsigned"))
            toks[iU] = " ".join((toks[iU], toks[iU+1]))
            del toks[iU+1]
        if len(toks) == 2:
            return toks[0].rstrip("&")
        else:
            nOpen = 0
            i = 0
            while i < len(toks) and ( i == 0 or nOpen != 0 ):
                nOpen += ( toks[i].count("<") - toks[i].count(">") )
                i += 1
            return " ".join(toks[:i]).rstrip("&")
    else:
        return "Float_t"

class CallMethod(TupleOp):
    """
    Call a method
    """
    __slots__ = ("name", "args", "_retType")
    def __init__(self, name, args, returnType=None, getFromRoot=True, canDefine=None):
        super(CallMethod, self).__init__()
        self.name = name ## NOTE can only be a hardcoded string this way
        self.args = tuple(adaptArg(arg) for arg in args)
        self._retType = returnType if returnType else CallMethod._initReturnType(name, getFromRoot=getFromRoot)
        self.canDefine = canDefine if canDefine is not None else all(a.canDefine for a in self.args)
    @staticmethod
    def _initReturnType(name, getFromRoot=True):
        mp = None
        if getFromRoot:
            try:
                from .root import gbl
                if "::" in name:
                    res = gbl
                    for tok in name.split("::"):
                        res = getattr(res, tok)
                    if res != gbl:
                        mp = res
                else:
                    mp = getattr(gbl, name)
            except Exception as ex:
                logger.error("Exception in getting method pointer {0}: {1}".format(name, ex))
        return guessReturnType(mp)
    def _clone(self, memo, select):
        argsCl = tuple(a.clone(memo=memo, select=select) for a in self.args)
        if select(self) or any(id(argCl) != id(arg) for argCl, arg in zip(argsCl, self.args)):
            return self.__class__(self.name, argsCl, returnType=self._retType, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in self.args:
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self._retType, self)
    def _repr(self):
        return "{0}({1!r}, ({2}), returnType={3!r})".format(self.__class__.__name__, self.name, ", ".join(repr(arg) for arg in self.args), self._retType)
    def _hash(self):
        return hash(tuple(chain([ self.__class__.__name__, self.name, self._retType], (hash(a) for a in self.args))))
    def _eq(self, other):
        return self.name == other.name and self._retType == other._retType and self.args == other.args
    # backends
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}({1})".format(self.name, ", ".join(defCache(arg) for arg in self.args))

class CallMemberMethod(TupleOp):
    """ Call a member method """
    __slots__ = ("this", "name", "args", "_retType")
    def __init__(self, this, name, args, returnType=None, canDefine=None):
        super(CallMemberMethod, self).__init__()
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        self.args = tuple(adaptArg(arg) for arg in args)
        self._retType = returnType if returnType else guessReturnType(getattr(this._typ, self.name))
        self.canDefine = canDefine if canDefine is not None else self.this.canDefine and all(a.canDefine for a in self.args)
    def _clone(self, memo, select):
        thisCl = self.this.clone(memo=memo, select=select)
        argsCl =  tuple(a.clone(memo=memo, select=select) for a in self.args)
        if select(self) or id(thisCl) != id(self.this) or any(id(argCl) != id(arg) for argCl, arg in zip(argsCl, self.args)):
            return self.__class__(thisCl, self.name, argsCl, returnType=self._retType, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in chain((self.this,), self.args):
                if select(arg):
                    yield arg
                yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self._retType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, ({3}), returnType={4!r})".format(self.__class__.__name__, self.this, self.name, ", ".join(repr(arg) for arg in self.args), self._retType)
    def _hash(self):
        return hash(tuple(chain([ self.__class__.__name__, self.this, self.name, self._retType ], (hash(a) for a in self.args))))
    def _eq(self, other):
        return self.this == other.this and self.name == other.name and self._retType == other._retType and self.args == other.args
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}({2})".format(defCache(self.this), self.name, ", ".join(defCache(arg) for arg in self.args))

class GetDataMember(TupleOp):
    """ Get a data member """
    __slots__ = ("this", "name")
    def __init__(self, this, name, canDefine=None):
        super(GetDataMember, self).__init__()
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        self.canDefine = canDefine if canDefine is not None else self.this.canDefine
    def _clone(self, memo, select):
        thisCl = self.this.clone(memo=memo, select=select)
        if select(self) or id(thisCl) != id(self.this):
            return self.__class__(thisCl, self.name, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            if select(self.this):
                yield self.this
            yield from self.this.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        if not self.name.startswith("_"):
            try:
                protoTp = self.this.result._typ
                proto = protoTp() ## should *in principle* work for most ROOT objects
                att = getattr(proto, self.name)
                tpNm = type(att).__name__
                if protoTp.__name__.startswith("pair<") and self.name in ("first", "second"):
                    tpNms = tuple(tok.strip() for tok in protoTp.__name__[5:-1].split(","))
                    return makeProxy((tpNms[0] if self.name == "first" else tpNms[1]), self)
                return makeProxy(tpNm, self)
            except Exception as e:
                print("Problem getting type of data member {0} of {1!r}".format(self.name, self.this), e)
        return makeProxy("void", self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.this, self.name)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.this), self.name))
    def _eq(self, other):
        return self.this == other.this and self.name == other.name
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}".format(defCache(self.this), self.name)

class ExtVar(TupleOp):
    """ Externally-defined variable (used by name) """
    __slots__ = ("typeName", "name")
    def __init__(self, typeName, name):
        super(ExtVar, self).__init__()
        self.typeName = typeName
        self.name = name
    def _clone(self, memo, select):
        if select(self):
            return self.__class__(self.typeName, self.name)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.typeName, self.name)
    def _hash(self):
        return hash(self._repr())
    def _eq(self, other):
        return self.typeName == other.typeName and self.name == other.name
    def get_cppStr(self, defCache=None):
        return self.name

class DefinedVar(TupleOp):
    """ Defined variable (used by name), first use will trigger definition """
    __slots__ = ("typeName", "definition", "_nameHint")
    def __init__(self, typeName, definition, nameHint=None):
        super(DefinedVar, self).__init__()
        self.typeName = typeName
        self.definition = definition
        self._nameHint = nameHint
    def _clone(self, memo, select):
        if select(self):
            return self.__class__(self.typeName, self.definition, nameHint=self._nameHint)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, nameHint={3!r})".format(self.__class__.__name__, self.typeName, self.definition, self._nameHint)
    def _hash(self):
        return hash((self.__class__.__name__, self.typeName, hash(self.definition), self._nameHint))
    def _eq(self, other):
        return self.typeName == other.typeName and self.definition == other.definition and self._nameHint == other._nameHint
    def get_cppStr(self, defCache=cppNoRedir):
        return defCache.symbol(self.definition, nameHint=self._nameHint)

class InitList(TupleOp):
    """ Initializer list """
    __slots__ = ("typeName", "elms")
    def __init__(self, typeName, elms, elmType=None, canDefine=None):
        super(InitList, self).__init__()
        self.typeName = typeName
        self.elms = tuple(adaptArg(e, typeHint=elmType) for e in elms)
        self.canDefine = canDefine if canDefine is not None else all(elm.canDefine for elm in self.elms)
    def _clone(self, memo, select):
        elmsCl = tuple(elm.clone(memo=memo, select=select) for elm in self.elms)
        if select(self) or any(id(elmCl) != id(elm) for elmCl, elm in zip(elmsCl, self.elms)):
            return self.__class__(self.typeName, elmsCl, canDefine=self.canDefine)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for elm in self.elms:
                if select(elm):
                    yield elm
                yield from elm.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def _repr(self):
        return "{0}<{1}>({2})".format(self.__class__.__name__, self.typeName, ", ".join(repr(elm) for elm in self.elms))
    def _hash(self):
        return hash((chain([ self.__class__.__name__, self.typeName], (hash(e) for e in self.elms))))
    def _eq(self, other):
        return self.typeName == other.typeName and self.elms == other.elms
    def get_cppStr(self, defCache=cppNoRedir):
        return "{{ {0} }}".format(", ".join(defCache(elm) for elm in self.elms))

class LocalVariablePlaceholder(TupleOp):
    """ Placeholder type for a local variable connected to an index (first step in a specific-to-general strategy) """
    __slots__ = ("typeHint", "_parent", "i")
    def __init__(self, typeHint, parent=None, i=None):
        super(LocalVariablePlaceholder, self).__init__()
        self.typeHint = typeHint
        self._parent = parent
        self.i = i ## FIXME this one is set **late** - watch out with what we call
        self.canDefine = False
    def _clone(self, memo, select):
        if select(self):
            return self.__class__(self.typeHint, parent=self._parent, i=self.i)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeHint, self)
    @property
    def name(self):
        if self.i is None:
            raise RuntimeError("Using LocalVariablePlaceholder before giving it an index")
        return "i{0:d}".format(self.i)
    def get_cppStr(self, defCache=None):
        return self.name
    def _repr(self):
        return "{0}({1!r}, i={2!r})".format(self.__class__.__name__, self.typeHint, self.i)
    def _hash(self):
        return hash((self.__class__.__name__, self.typeHint, self.i))
    def _eq(self, other):
        ## NOTE this breaks the infinite recursion, but may not be 100% safe
        ## what should save the nested case is that the repr(parent) will be different for different levels of nesting
        ## since all LVP's are supposed to have an index, confusion between cases where they are combined differently should be eliminated as well
        return self.typeHint == other.typeHint and self.i == other.i and ( id(self._parent) == id(other._parent) or ( hash(self._parent) == hash(other._parent) and self._parent.__class__ == other._parent.__class__ and repr(self._parent) == repr(other._parent) ) )

def collectNodes(expr, defCache=cppNoRedir, select=(lambda nd : True), stop=(lambda nd : False), includeLocal=True):
    # simple helper
    if select(expr):
        yield expr
    if not stop(expr):
        yield from expr.deps(select=select, includeLocal=includeLocal)

def _collectDeps(self, defCache=cppNoRedir):
    ## first pass (will trigger definitions, if necessary)
    for dep in self.deps(defCache=defCache, select=lambda op : defCache.shouldDefine(op)):
        cn = defCache(dep)
        if not cn:
            logger.warning("Probably a problem in triggering definition for {0}".format(dep))
    return set(dep for dep in self.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn) or isinstance(op, LocalVariablePlaceholder) or defCache.shouldDefine(op) or (defCache._getColName(op) is not None))))

def _canDefine(expr, local):
    return not any(ind is not None for ind in collectNodes(expr, select=(lambda nd : isinstance(nd, LocalVariablePlaceholder) and nd not in local), stop=(lambda nd : nd.canDefine), includeLocal=False))

def _convertFunArgs(deps, defCache=cppNoRedir):
    capDeclCall = []
    for ld in deps:
        if isinstance(ld, GetArrayColumn):
            capDeclCall.append((
                "&{0}".format(ld.name),
                "const ROOT::VecOps::RVec<{0}>& {1}".format(ld.typeName, ld.name),
                ld.name))
        elif isinstance(ld, GetColumn):
            capDeclCall.append((
                "&{0}".format(ld.name),
                "const {0}& {1}".format(ld.typeName, ld.name),
                ld.name))
        elif isinstance(ld, LocalVariablePlaceholder):
            if not ld.name:
                print("ERROR: no name for local {0}".format(ld))
            capDeclCall.append((
                ld.name,
                "{0} {1}".format(ld.typeHint, ld.name),
                ld.name))
        elif defCache.shouldDefine(ld) or (defCache._getColName(ld) is not None):
            nm = defCache._getColName(ld)
            if not nm:
                print("ERROR: no column name for {0}".format(ld))
            if not any("&{0}".format(nm) == icap for icap,idecl,icall in capDeclCall):
                capDeclCall.append((
                    "&{0}".format(nm),
                    "const {0}& {1}".format(ld.result._typeName, nm),
                    nm))
            else:
                print("WARNING: dependency {0} is there twice".format(nm))
        else:
            raise AssertionError("Dependency with unknown type: {0}".format(ld))
    return zip(*sorted(capDeclCall, key=(lambda elm : elm[1]))) ## sort by declaration (alphabetic for the same type)

def _normFunArgs(expr, args, argNames):
    newExpr = expr
    newArgs = args
    for i,argN in sorted(list(enumerate(argNames)), key=(lambda elm : len(elm[1])), reverse=True):
        newName = "myArg{0:d}".format(i, argN)
        if sum(1 for ia in newArgs if argN in ia) != 1: ## should be in one and only one argument
            raise RuntimeError("{0} is in more than one (or none) of {1}".format(argN, newArgs))
        newArgs = [ (ia.replace(argN, newName) if argN in ia else ia) for ia in newArgs ]
        newExpr = newExpr.replace(argN, newName)
    return newExpr, newArgs

class Select(TupleOp):
    """ Define a selection on a range """
    __slots__ = ("rng", "predExpr", "_i")
    def __init__(self, rng, predExpr, idx, canDefine=None):
        super(Select, self).__init__()
        self.rng = rng ## input indices
        self.predExpr = predExpr
        self._i = idx
        self.canDefine = canDefine if canDefine is not None else self.rng.canDefine and _canDefine(self.predExpr, (self._i,))
    def _clone(self, memo, select):
        rngCl = self.rng.clone(memo=memo, select=select)
        predCl = self.predExpr.clone(memo=memo, select=select)
        iCl = self._i.clone(memo=memo, select=select)
        if select(self) or id(rngCl) != id(self.rng) or id(predCl) != id(self.predExpr) or id(iCl) != id(self._i):
            return self.__class__(rngCl, predCl, iCl, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(rng, pred):
        """ Factory method from a range and predicate (callable) """
        idx = LocalVariablePlaceholder(SizeType)
        predExpr = adaptArg(pred(rng._getItem(idx.result)))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(predExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Select(adaptArg(rng._idxs), predExpr, idx)
        idx._parent = res
        from .treeproxies import SelectionProxy
        return SelectionProxy(rng._base, res, valueType=rng.valueType)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.predExpr):
                if select(arg) and ( includeLocal or arg != self._i ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(SizeType), itemType=SizeType)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.predExpr, self._i)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.rng), hash(self.predExpr), hash(self._i)))
    def _eq(self, other):
        return self.rng == other.rng and self.predExpr == other.predExpr and self._i == other._i
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::select({idxs},\n    [{captures}] ( {i} ) {{ return {predExpr}; }})".format(
                idxs=defCache(self.rng),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                predExpr=defCache(self.predExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType="ROOT::VecOps::RVec<{0}>".format(SizeType), args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Sort(TupleOp):
    """ Sort a range (ascendingly) by the value of a function on each element """
    __slots__ = ("rng", "funExpr", "_i")
    def __init__(self, rng, funExpr, idx, canDefine=None):
        super(Sort, self).__init__()
        self.rng = rng ## input indices
        self.funExpr = funExpr
        self._i = idx
        self.canDefine = canDefine if canDefine is not None else self.rng.canDefine and _canDefine(self.funExpr, (self._i,))
    def _clone(self, memo, select):
        rngCl = self.rng.clone(memo=memo, select=select)
        funCl = self.funExpr.clone(memo=memo, select=select)
        iCl = self._i.clone(memo=memo, select=select)
        if select(self) or id(rngCl) != id(self.rng) or id(funCl) != id(self.funExpr) or id(iCl) != id(self._i):
            return self.__class__(rngCl, funCl, iCl, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(rng, fun):
        idx = LocalVariablePlaceholder(SizeType)
        funExpr = adaptArg(fun(rng._getItem(idx.result)))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(funExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Sort(adaptArg(rng._idxs), funExpr, idx)
        idx._parent = res
        from .treeproxies import SelectionProxy
        return SelectionProxy(rng._base, res, valueType=rng.valueType)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.funExpr):
                if select(arg) and ( includeLocal or arg != self._i ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(SizeType), itemType=SizeType)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.funExpr, self._i)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.rng), hash(self.funExpr), hash(self._i)))
    def _eq(self, other):
        return self.rng == other.rng and self.funExpr == other.funExpr and self._i == other._i
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::sort({idxs},\n    [{captures}] ( {i} ) {{ return {funExpr}; }})".format(
                idxs=defCache(self.rng),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                funExpr=defCache(self.funExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType="ROOT::VecOps::RVec<{0}>".format(SizeType), args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Map(TupleOp):
    """ Create a list of derived values for a collection (mostly useful for storing on skims) """
    __slots__ = ("rng", "funExpr", "_i", "typeName")
    def __init__(self, rng, funExpr, idx, typeName, canDefine=None):
        super(Map, self).__init__()
        self.rng = rng ## input indices
        self.funExpr = funExpr
        self._i = idx
        self.typeName = typeName
        self.canDefine = canDefine if canDefine is not None else self.rng.canDefine and _canDefine(self.funExpr, (self._i,))
    def _clone(self, memo, select):
        rngCl = self.rng.clone(memo=memo, select=select)
        funCl = self.funExpr.clone(memo=memo, select=select)
        iCl = self._i.clone(memo=memo, select=select)
        if select(self) or id(rngCl) != id(self.rng) or id(funCl) != id(self.funExpr) or id(iCl) != id(self._i):
            return self.__class__(rngCl, funCl, iCl, self.typeName, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(rng, fun, typeName=None):
        idx = LocalVariablePlaceholder(SizeType)
        val = fun(rng._getItem(idx.result))
        funExpr = adaptArg(val)
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(funExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Map(adaptArg(rng._idxs), funExpr, idx, typeName=(typeName if typeName is not None else val._typeName))
        idx._parent = res
        return res.result
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.funExpr):
                if select(arg) and ( includeLocal or arg != self._i ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, typeName="ROOT::VecOps::RVec<{0}>".format(self.typeName), itemType=self.typeName)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self.rng, self.funExpr, self._i, self.typeName)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.rng), hash(self.funExpr), hash(self._i), self.typeName))
    def _eq(self, other):
        return self.rng == other.rng and self.funExpr == other.funExpr and self._i == other._i and self.typeName == other.typeName
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::map<{valueType}>({idxs},\n    [{captures}] ( {i} ) {{ return {funExpr}; }})".format(
                valueType=self.typeName,
                idxs=defCache(self.rng),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                funExpr=defCache(self.funExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType="ROOT::VecOps::RVec<{0}>".format(self.typeName), args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Next(TupleOp):
    """ Define a search (first matching item, for a version that processes the whole range see Reduce) """
    __slots__ = ("rng", "predExpr", "_i")
    def __init__(self, rng, predExpr, idx, canDefine=None):
        super(Next, self).__init__()
        self.rng = rng ## input indices
        self.predExpr = predExpr
        self._i = idx
        self.canDefine = canDefine if canDefine is not None else self.rng.canDefine and _canDefine(self.predExpr, (self._i,))
    def _clone(self, memo, select):
        rngCl = self.rng.clone(memo=memo, select=select)
        predCl = self.predExpr.clone(memo=memo, select=select)
        iCl = self._i.clone(memo=memo, select=select)
        if select(self) or id(rngCl) != id(self.rng) or id(predCl) != id(self.predExpr) or id(iCl) != id(self._i):
            return self.__class__(rngCl, predCl, iCl, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(rng, pred): ## FIXME you are here
        idx = LocalVariablePlaceholder(SizeType)
        predExpr = adaptArg(pred(rng._getItem(idx.result)))
        idx.i = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(predExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))+1
        res = Next(adaptArg(rng._idxs), predExpr, idx)
        idx._parent = res
        return rng._getItem(res)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.predExpr):
                if select(arg) and ( includeLocal or arg != self._i ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(SizeType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.rng, self.predExpr, self._i)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.rng), hash(self.predExpr), hash(self._i)))
    def _eq(self, other):
        return self.rng == other.rng and self.predExpr == other.predExpr and self._i == other._i
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::next({idxs},\n     [{captures}] ( {i} ) {{ return {predexpr}; }}, -1)".format(
                idxs=defCache(self.rng),
                captures=", ".join(captures),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                predexpr=defCache(self.predExpr),
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType=SizeType, args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Reduce(TupleOp):
    """ Reduce a range to a value (could be a transformation, index...) """
    def __init__(self, rng, resultType, start, accuExpr, idx, prevRes, canDefine=None):
        super(Reduce, self).__init__()
        self.rng = rng ## input indices
        self.resultType = resultType
        self.start = start
        self.accuExpr = accuExpr
        self._i = idx
        self._prevRes = prevRes
        self.canDefine = canDefine if canDefine is not None else self.rng.canDefine and _canDefine(start, (self._i, self._prevRes)) and _canDefine(accuExpr, (self._i, self._prevRes))
    def _clone(self, memo, select):
        rngCl = self.rng.clone(memo=memo, select=select)
        startCl = self.start.clone(memo=memo, select=select)
        accuCl = self.accuExpr.clone(memo=memo, select=select)
        iCl = self._i.clone(memo=memo, select=select)
        prevCl = self._prevRes.clone(memo=memo, select=select)
        if select(self) or id(rngCl) != id(self.rng) or id(startCl) != id(self.start) or id(accuCl) != id(self.accuExpr) or id(iCl) != id(self._i) or id(prevCl) != id(self._prevRes):
            return self.__class__(rngCl, self.resultType, startCl, accuCl, iCl, prevCl, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(rng, start, accuFun):
        resultType = start._typeName
        idx = LocalVariablePlaceholder(SizeType)
        prevRes = LocalVariablePlaceholder(resultType, i=-1)
        accuExpr = adaptArg(accuFun(prevRes.result, rng._getItem(idx.result)))
        maxLVIdx = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(accuExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))
        idx.i = maxLVIdx+1
        prevRes.i = maxLVIdx+2

        res = Reduce(adaptArg(rng._idxs), resultType, adaptArg(start), accuExpr, idx, prevRes)
        idx._parent = res
        prevRes._parent = res
        return res.result
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.start, self.accuExpr):
                if select(arg) and ( includeLocal or arg not in (self._i, self._prevRes) ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in (self._i, self._prevRes):
                        yield dp
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.resultType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r}, {5!r}, {6!r})".format(self.__class__.__name__, self.rng, self.resultType, self.start, self.accuExpr, self._i, self._prevRes)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.rng), self.resultType, hash(self.start), hash(self.accuExpr), hash(self._i), hash(self._prevRes)))
    def _eq(self, other):
        return self.rng == other.rng and self.resultType == other.resultType and self.start == other.start and self.accuExpr == other.accuExpr and self._i == other._i and self._prevRes == other._prevRes
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = "rdfhelpers::reduce({idxs}, {start},\n     [{captures}] ( {prevRes}, {i} ) {{ return {accuexpr}; }})".format(
                idxs=defCache(self.rng),
                start=defCache(self.start),
                captures=", ".join(captures),
                prevRes="{0} {1}".format(self._prevRes.typeHint, self._prevRes.name),
                i="{0} {1}".format(self._i.typeHint, self._i.name),
                accuexpr=defCache(self.accuExpr)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType=self.resultType, args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

class Combine(TupleOp):
    __slots__ = ("ranges", "candPredExpr", "_i")
    def __init__(self, ranges, candPredExpr, idx, canDefine=None):
        super(Combine, self).__init__()
        self.ranges = ranges ## [ input index lists ]
        self.candPredExpr = candPredExpr
        self._i = idx
        self.canDefine = canDefine if canDefine is not None else all(rng.canDefine for rng in self.ranges) and _canDefine(candPredExpr, self._i)
    @property
    def n(self):
        return len(self.ranges)
    def _clone(self, memo, select):
        rngCl = tuple(rng.clone(memo=memo, select=select) for rng in self.ranges)
        predCl = self.candPredExpr.clone(memo=memo, select=select)
        iCl = tuple(i.clone(memo=memo, select=select) for i in self._i)
        if select(self) or any(id(rCl) != id(rng) for rCl, rng in zip(rngCl, self.ranges)) or id(predCl) != id(self.candPredExpr) or any(id(iiCl) != id(iOr) for iiCl,iOr in zip(iCl, self._i)):
            return self.__class__(rngCl, predCl, iCl, canDefine=self.canDefine)
    @staticmethod
    def fromRngFun(num, ranges, candPredFun, samePred=None):
        ranges = ranges if len(ranges) > 1 else list(repeat(ranges[0], num))
        idx = tuple(LocalVariablePlaceholder(SizeType, i=-1-i) for i in range(num))
        candPred = candPredFun(*( rng._getItem(iidx.result) for rng,iidx in zip(ranges, idx)))
        candPredExpr = adaptArg(candPred)
        if samePred:
            from . import treefunctions as op
            areDiff = op.AND(*(samePred(ra._getItem(ia.result), rb._getItem(ib.result))
                    for ((ia, ra), (ib, rb)) in combinations(zip(idx, ranges), 2)
                    if ra._base == rb._base))
            if len(areDiff.op.args) > 0:
                candPredExpr = adaptArg(op.AND(areDiff, candPred))
        maxLVIdx = max(chain([-1], ((nd.i if nd.i is not None else -1) for nd in collectNodes(candPredExpr,
            select=(lambda nd : isinstance(nd, LocalVariablePlaceholder))))))
        for i,ilvp in enumerate(idx):
            ilvp.i = maxLVIdx+1+i
        res = Combine(tuple(adaptArg(rng._idxs) for rng in ranges), candPredExpr, idx)
        for ilvp in idx:
            ilvp._parent = res
        from .treeproxies import CombinationListProxy
        return CombinationListProxy(ranges, res)
    @property
    def resultType(self):
        return "ROOT::VecOps::RVec<rdfhelpers::Combination<{0:d}>>".format(self.n)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in chain(self.ranges, [self.candPredExpr]):
                if select(arg) and ( includeLocal or arg not in self._i ):
                    yield arg
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.resultType, self)
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.ranges, self.candPredExpr, self._i)
    def _hash(self):
        return hash((self.__class__.__name__, tuple(hash(r) for r in self.ranges), hash(self.candPredExpr), tuple(hash(i) for i in self._i)))
    def _eq(self, other):
        return self.ranges == other.ranges and self.candPredExpr == other.candPredExpr and self._i == other._i
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(self, defCache=defCache)
        captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
        expr = ("rdfhelpers::combine(\n"
            "     [{captures}] ( {predIdxArgs} ) {{ return {predExpr}; }},\n"
            "     {ranges})").format(
                captures=", ".join(captures),
                predIdxArgs=", ".join("{0} {1}".format(i.typeHint, i.name) for i in self._i),
                predExpr = defCache(self.candPredExpr),
                ranges=", ".join(defCache(rng) for rng in self.ranges)
                )
        if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
            return expr
        else:
            expr_n, paramDecl_n = _normFunArgs(expr, paramDecl, paramCall)
            funName = defCache.symbol(expr_n, resultType=self.resultType, args=", ".join(paramDecl_n))
            return "{0}({1})".format(funName, ", ".join(paramCall))

## FIXME to be implemented
class PseudoRandom(TupleOp):
    """ Pseudorandom number (integer or float) within range """
    def __init__(self, xMin, xMax, seed, isIntegral=False):
        super(PseudoRandom, self).__init__()
        self.xMin = xMin
        self.xMax = xMax
        self.seed = seed
        self.isIntegral = isIntegral
    @property
    def resultType(self):
        return "Int_" if self.isIntegral else "Float_t"
    ## deps from xMin, xMax and seed
    ## seed can be event-based or object-based, depending?
    ## TODO implement C++ side as well

class OpWithSyst(ForwardingOp):
    """ Interface and base class for nodes that can change the systematic variation of something they wrap """
    def __init__(self, wrapped, systName, variations=None):
        super(OpWithSyst, self).__init__(wrapped)
        self.systName = systName
        self.variations = variations
    def _clone(self, memo, select):
        clWr = self.wrapped.clone(memo=memo, select=select)
        if select(self) or id(clWr) != id(self.wrapped):
            return self.__class__(clWr, self.systName, variations=self.variations)
    def changeVariation(self, newVariation):
        """ Assumed to be called on a fresh copy - *will* change the underlying value
        """
        if self._cache: # validate this assumption
            raise RuntimeError("Cannot change variation of an expression that is already frozen")
        if newVariation not in self.variations:
            raise ValueError("Invalid variation: {0}".format(newVariation))
        self._changeVariation(newVariation)
    def _changeVariation(self, newVariation):
        """ changeVariation specifics (after validating newVariation, and changability) - to be implemented by concrete classes """
        pass
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self.wrapped, self.systName, self.variations)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.wrapped), self.systName, self.variations))
    def _eq(self, other):
        return super(OpWithSyst, self)._eq(other) and self.systName == other.systName and self.variations == other.variations

class ScaleFactorWithSystOp(OpWithSyst):
    """ Scalefactor (ILeptonScaleFactor::get() call), to be modified with Up/Down variations (these are cached) """
    def __init__(self, wrapped, systName, variations=None):
        super(ScaleFactorWithSystOp, self).__init__(wrapped, systName, variations=(variations if variations else tuple("{0}{1}".format(systName, vard) for vard in ["up", "down"])))
    def _changeVariation(self, newVariation):
        """ Assumed to be called on a fresh copy - *will* change the underlying value """
        newVariation = (newVariation[len(self.systName):] if newVariation.startswith(self.systName) else newVariation).capitalize() ## translate to name in C++
        wrOrig = self.wrapped
        assert isinstance(wrOrig, CallMemberMethod) and isinstance(wrOrig.args[-1], ExtVar)
        if wrOrig.args[-1].name == "Nominal" and newVariation != "Nominal":
            replArgs = tuple(list(wrOrig.args)[:-1]+[ ExtVar(wrOrig.args[-1].typeName, newVariation) ])
            self.wrapped = CallMemberMethod(wrOrig.this, wrOrig.name, replArgs, returnType=wrOrig._retType, canDefine=wrOrig.canDefine)

class SystAltOp(OpWithSyst):
    """ Change the wrapped operation (from a map) """
    def __init__(self, wrapped, name, varMap, valid=None):
        super(SystAltOp, self).__init__(wrapped, name)
        self.variations = valid if valid is not None else tuple(varMap.keys())
        self.varMap = varMap
    def _clone(self, memo, select):
        clWr = self.wrapped.clone(memo=memo, select=select)
        varCl = dict((nm, vop.clone(memo=memo, select=select)) for nm,vop in self.varMap.items())
        if select(self) or id(clWr) != id(self.wrapped) or any(id(vCl) != id(self.varMap[ky]) for ky,vCl in varCl.items()):
            return self.__class__(clWr, self.systName, varCl, valid=tuple(self.variations))
    def _repr(self):
        return "{0}({1!r}, {2!r}, {3!r}, {4!r})".format(self.__class__.__name__, self.wrapped, self.systName, self.variations, self.varMap)
    def _hash(self):
        return hash((self.__class__.__name__, hash(self.wrapped), self.systName, tuple(self.variations), tuple((ky, hash(val)) for ky,val in self.varMap.items())))
    def _eq(self, other):
        return super(SystAltOp, self)._eq(other) and self.variations == other.variations and self.varMap == other.varMap
    def _changeVariation(self, newVariation):
        if newVariation in self.varMap:
            self.wrapped = self.varMap[newVariation]

def collectSystVars(exprs):
    systVars = {}
    for sfs in chain.from_iterable(collectNodes(expr, select=(lambda nd : isinstance(nd, OpWithSyst) and nd.systName and nd.variations)) for expr in exprs):
        if sfs.systName not in systVars:
            systVars[sfs.systName] = list(sfs.variations)
        else:
            for sv in sfs.variations:
                if sv not in systVars[sfs.systName]:
                    systVars[sfs.systName].append(sv)
    return systVars
def mergeSystVars(svA, svB):
    ## returns A updated with B
    for systName, variations in svB.items():
        if systName not in svA:
            svA[systName] = list(variations)
        else:
            for sv in variations:
                if sv not in svA[systName]:
                    svA[systName].append(sv)
    return svA
