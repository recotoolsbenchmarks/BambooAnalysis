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

class TupleOp(object):
    """ Interface & base class for operations on leafs and resulting objects / values """
    def deps(self, defCache=None, select=(lambda x : True), includeLocal=False):
        yield from []
    @property
    def result(self):
        pass
    ## subclasses should define __eq__, __repr__ and __hash__

## implementations are split out, see treeproxies
class TupleBaseProxy(object):
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

class CppStrRedir(object):
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

cppNoRedir = CppStrRedir()

SizeType = "std::size_t"

class Const(TupleOp):
    """ Hard-coded number (or expression) """
    __slots__ = ("typeName", "value")
    def __init__(self, typeName, value):
        self.typeName = typeName
        self.value = value
        super(Const, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __repr__(self):
        return "Const({0!r})".format(self.value)
    def __eq__(self, other):
        return isinstance(other, Const) and ( self.typeName == other.typeName ) and ( self.value == other.value )
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
    def __init__(self, typeName, name):
        self.typeName = typeName
        self.name = name
        super(GetColumn, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, GetColumn) and self.name == other.name and self.typeName == other.typeName
    def __repr__(self):
        return "GetColumn({0!r}, {1!r})".format(self.typeName, self.name)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=None):
        return self.name

class GetArrayColumn(TupleOp):
    """ Get the number from a leaf """
    def __init__(self, typeName, name, length):
        self.typeName = typeName
        self.name = name
        self.length = length
        super(GetArrayColumn, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if select(self.length):
            yield self.length
        yield from self.length.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self, makeProxy(SizeType, self.length))
    def __eq__(self, other):
        return isinstance(other, GetArrayColumn) and ( self.typeName == other.typeName ) and ( self.name == other.name ) and ( self.length == other.length )
    def __repr__(self):
        return "GetArrayColumn({0!r}, {1!r}, {2!r})".format(self.typeName, self.name, self.length)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=None):
        return self.name

## helper
def adaptArg(arg, typeHint=None):
    from .treeproxies import TupleBaseProxy
    if isinstance(arg, TupleBaseProxy):
        return arg.op
    elif isinstance(arg, TupleOp):
        return arg
    elif typeHint is not None:
        if str(arg) == arg: ## string, needs quote
            return Const(typeHint, '"{}"'.format(arg))
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
    #
    , "abs"   : lambda cppStr,arg : "std::abs( {0} )".format(cppStr(arg))
    , "sqrt"  : lambda cppStr,arg : "std::sqrt( {0} )".format(cppStr(arg))
    , "pow"   : lambda cppStr,a1,a2 : "std::pow( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "exp"   : lambda cppStr,arg : "std::exp( {0} )".format(cppStr(arg))
    , "log"   : lambda cppStr,arg : "std::log( {0} )".format(cppStr(arg))
    , "log10" : lambda cppStr,arg : "std::log10( {0} )".format(cppStr(arg))
    , "max"   : lambda cppStr,a1,a2 : "std::max( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    , "min"   : lambda cppStr,a1,a2 : "std::min( {0}, {1} )".format(cppStr(a1), cppStr(a2))
    #
    , "switch" : lambda cppStr,test,trueBr,falseBr : "( {0} ) ? ( {1} ) : ( {2} )".format(cppStr(test), cppStr(trueBr), cppStr(falseBr))
    }

class MathOp(TupleOp):
    """ Mathematical function N->1, e.g. sin, abs, ( lambda x, y : x*y ) """
    def __init__(self, op, *args, **kwargs):
        self.outType = kwargs.pop("outType", "Double_t")
        assert len(kwargs) == 0
        self.op = op
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
        super(MathOp, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.outType, self)
    def __eq__(self, other):
        return isinstance(other, MathOp) and ( self.outType == other.outType ) and ( self.op == other.op ) and ( len(self.args) == len(other.args) ) and all( ( sa == oa ) for sa, oa in zip(self.args, other.args))
    def __repr__(self):
        return "MathOp({0}, {1})".format(self.op, ", ".join(repr(arg) for arg in self.args))
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return mathOpFuns_cppStr[self.op](defCache, *self.args)

class GetItem(TupleOp):
    """ Get item from array (from function call or from array leaf) """
    def __init__(self, arg, valueType, index, indexType=SizeType):
        self.arg = adaptArg(arg)
        self.typeName = valueType
        self._index = adaptArg(index, typeHint=SizeType)
        super(GetItem, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
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
    def __eq__(self, other):
        return isinstance(other, GetItem) and ( self.arg == other.arg ) and ( self.typeName == other.typeName ) and ( self._index == other._index )
    def __repr__(self):
        return "GetItem({0!r}, {1!r})".format(self.arg, self._index)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}[{1}]".format(defCache(self.arg), defCache(self._index))

class Construct(TupleOp):
    def __init__(self, typeName, args):
        self.typeName = typeName
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
        super(Construct, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, Construct) and ( self.typeName == other.typeName ) and len(self.args) == len(other.args) and all( ( aa == ab ) for aa,ab in zip(self.args, other.args) )
    def __repr__(self):
        return "Construct({0!r}, {1})".format(self.typeName, ", ".join(repr(a) for a in self.args))
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}{{{1}}}".format(self.typeName, ", ".join(defCache(a) for a in self.args))

def guessReturnType(mp):
    toks = list(mp.func_doc.split())
    ## left and right strip const * and &
    while toks[-1].rstrip("&") in ("", "const"):
        toks = toks[:-1]
    while toks[0].rstrip("&") in ("", "const"):
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

class CallMethod(TupleOp):
    """
    Call a method
    """
    def __init__(self, name, args):
        self.name = name ## NOTE can only be a hardcoded string this way
        self._mp = None
        try:
            from cppyy import gbl
            if "::" in name:
                res = gbl
                for tok in name.split("::"):
                    res = getattr(res, tok)
                if res != gbl:
                    self._mp = res
            else:
                self._mp  = getattr(gbl, name)
        except Exception as ex:
            logger.error("Exception in getting method pointer: {0}".format(ex))
            self._mp = None
        self.args = tuple(adaptArg(arg) for arg in args)
        super(CallMethod, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        retTypeN = "Float_t"
        if self._mp and hasattr(self._mp, "func_doc") and hasattr(self._mp, "func_name"):
            fdoc = self._mp.func_doc
            fnam = self._mp.func_name
            toks = fdoc.split()
            if toks[0] == "const":
                del toks[0]
            out = ""
            for tok in toks:
                if fnam in tok:
                    break
                out += " {0}".format(tok)
            retTypeN = out.strip().strip("*&")
        from .treeproxies import makeProxy
        return makeProxy(retTypeN, self)
    def __eq__(self, other):
        return isinstance(other, CallMethod) and ( self.name == other.name ) and ( len(self.args) == len(other.args) ) and all( ( sa == oa ) for sa, oa in zip(self.args, other.args))
    def __repr__(self):
        return "CallMethod({0!r}, ({1}))".format(self.name, ", ".join(repr(arg) for arg in self.args))
    def __hash__(self):
        return hash(self.__repr__())
    # backends
    def get_cppStr(self, defCache=cppNoRedir):
        if not defCache.backend.shouldDefine(self, defCache=defCache):
            return "{0}({1})".format(self.name, ", ".join(defCache(arg) for arg in self.args))
        else: ## go through a symbol
            depList = _collectDeps(self.args, [], defCache=defCache)
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = "{name}({args})\n".format(name=self.name, args=", ".join(defCache(arg) for arg in self.args))
            funName = defCache.symbol(expr, resultType=self.result._typeName, args=paramDecl)
            return "{0}({1})".format(funName, paramCall)

class CallMemberMethod(TupleOp):
    """ Call a member method """
    def __init__(self, this, name, args):
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        self._mp  = getattr(this._typ, name)
        self.args = tuple(adaptArg(arg) for arg in args)
        super(CallMemberMethod, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for arg in chain((self.this,), self.args):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        retTypeN = guessReturnType(self._mp)
        from .treeproxies import makeProxy
        return makeProxy(retTypeN, self)
    def __eq__(self, other):
        return isinstance(other, CallMemberMethod) and ( self.this == other.this ) and ( self.name == other.name ) and ( len(self.args) == len(other.args) ) and all( ( sa == oa ) for sa, oa in zip(self.args, other.args))
    def __repr__(self):
        return "CallMemberMethod({0!r}, {1!r}, ({2}))".format(self.this, self.name, ", ".join(repr(arg) for arg in self.args))
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}({2})".format(defCache(self.this), self.name, ", ".join(defCache(arg) for arg in self.args))

class GetDataMember(TupleOp):
    """ Get a data member """
    __slots__ = ("this", "name")
    def __init__(self, this, name):
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        super(GetDataMember, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
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
    def __eq__(self, other):
        return isinstance(other, GetDataMember) and ( self.this == other.this ) and ( self.name == other.name )
    def __repr__(self):
        return "GetDataMember({0!r}, {1!r})".format(self.this, self.name)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}.{1}".format(defCache(self.this), self.name)

class ExtVar(TupleOp):
    """ Externally-defined variable (used by name) """
    def __init__(self, typeName, name):
        self.typeName = typeName
        self.name = name
        super(ExtVar, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, ExtVar) and ( self.typeName == other.typeName ) and ( self.name == other.name )
    def __repr__(self):
        return "ExtVar({0!r}, {1!r})".format(self.typeName, self.name)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=None):
        return self.name

class DefinedVar(TupleOp):
    """ Defined variable (used by name), first use will trigger definition """
    def __init__(self, typeName, definition, nameHint=None):
        self.typeName = typeName
        self.definition = definition
        self._nameHint = nameHint
        super(DefinedVar, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, DefinedVar) and ( self.typeName == other.typeName ) and ( self.definition == other.definition )
    def __repr__(self):
        return "DefinedVar({0!r}, {1!r})".format(self.typeName, self.definition)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return defCache.symbol(self.definition, nameHint=self._nameHint)

class InitList(TupleOp):
    def __init__(self, typeName, elmType, elms):
        self.typeName = typeName
        self.elms = tuple(adaptArg(e, typeHint=elmType) for e in elms)
        super(InitList, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        for elm in self.elms:
            if select(elm):
                yield elm
            yield from elm.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, InitList) and ( self.typeName == other.typeName ) and len(self.elms) == len(other.elms) and all( ( ea == eb ) for ea, eb in zip(self.elms, other.elms) )
    def __repr__(self):
        return "InitList<{0}>({1})".format(self.typeName, ", ".join(repr(elm) for elm in self.elms))
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{{ {0} }}".format(", ".join(defCache(elm) for elm in self.elms))

class LocalVariablePlaceholder(TupleOp):
    """ Placeholder type for a local variable connected to an index (first step in a specific-to-general strategy) """
    def __init__(self, typeHint, name=None):
        self.name = name
        self.typeHint = typeHint
        super(LocalVariablePlaceholder, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeHint, self)
    ## NOTE no __eq__, __repr__ and __hash__ because same name and type does not mean they are equal (should use reference)
    def get_cppStr(self, defCache=None):
        if not self.name:
            raise RuntimeError("Using LocalVariablePlaceholder before giving it a name")
        return self.name

def _collectDeps(exprs, ownLocal, defCache=cppNoRedir):
    ## first pass (will trigger definitions, if necessary)
    exprs1, exprs2 = tee(exprs, 2)
    for dep in chain.from_iterable(expr.deps(defCache=defCache, select=lambda op : defCache.backend.shouldDefine(op, defCache=defCache)) for expr in exprs1):
        cn = defCache(dep)
        if not cn:
            print("WARNING: Probably a problem in triggering definition for {0}".format(dep))
    return set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn)
                or defCache.backend.shouldDefine(op, defCache=defCache) or ( isinstance(op, LocalVariablePlaceholder) and op not in ownLocal )
                ))
            for expr in exprs2))

def _collectLocalVars(exprs, ownLocal, defCache=cppNoRedir):
    localVarsToName = set(chain.from_iterable(
        expr.deps(defCache=defCache, select=(lambda op : isinstance(op, LocalVariablePlaceholder) and not op.name), includeLocal=True)
        for expr in exprs))
    for lv in ownLocal:
        if not lv.name:
            localVarsToName.add(lv)
    return localVarsToName

@contextmanager
def _nameLocalVars(localVars):
    for i,lv in zip(count(), localVars):
        if lv.name:
            raise RuntimeError("Should name {0!r} {1}?".format(lv, "i{0:d".format(i)))
        lv.name = "i{0:d}".format(i)
    yield
    for lv in localVars:
        lv.name = None

def _convertFunArgs(deps, defCache=cppNoRedir):
    captures, paramDecl, paramCall = [], [], []
    for ld in deps:
        if isinstance(ld, GetArrayColumn):
            captures.append("&{0}".format(ld.name))
            paramDecl.append("const ROOT::VecOps::RVec<{0}>& {1}".format(ld.typeName, ld.name))
            paramCall.append(ld.name)
        elif isinstance(ld, GetColumn):
            captures.append("&{0}".format(ld.name))
            paramDecl.append("const {0}& {1}".format(ld.typeName, ld.name))
            paramCall.append(ld.name)
        elif isinstance(ld, LocalVariablePlaceholder):
            if not ld.name:
                print("ERROR: no name for local {0}".format(ld))
            captures.append(ld.name)
            paramDecl.append("{0} {1}".format(ld.typeHint, ld.name))
            paramCall.append(ld.name)
        elif defCache.backend.shouldDefine(ld, defCache=defCache):
            nm = defCache._getColName(ld)
            if not nm:
                print("ERROR: no column name for {0}".format(ld))
            if "&{0}".format(nm) not in captures:
                captures.append("&{0}".format(nm))
                paramDecl.append("const {0}& {1}".format(ld.result._typeName, nm))
                paramCall.append(nm)
            else:
                print("WARNING: dependency {0} is there twice".format(nm))
        else:
            raise AssertionError("Dependency with unknown type: {1}, or {0!r}".format(ld, cppStr))
    return ",".join(captures), ", ".join(paramDecl), ", ".join(paramCall)

class Select(TupleOp):
    """ Define a selection on a range """
    def __init__(self, rng, pred):
        self.rng = rng ## PROXY
        self._i = LocalVariablePlaceholder(SizeType)
        self.predExpr = adaptArg(pred(self.rng._base[self._i.result]))
        super(Select, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (adaptArg(self.rng), self.predExpr):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import VectorProxy
        return VectorProxy(self, "ROOT::VecOps::RVec<{0}>".format(SizeType))
    def __eq__(self, other):
        return isinstance(other, Select) and ( self.rng == other.rng ) and ( self.predExpr == other.predExpr )
    def __repr__(self):
        return "Select({0!r}, {1!r})".format(self.rng, self.predExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.predExpr), (self._i,), defCache=defCache)
        ## should only be a non-empty list for the topmost op in case of nested expressions
        localVarsToName = _collectLocalVars((self.rng, self.predExpr), (self._i,), defCache=defCache)
        with _nameLocalVars(localVarsToName):
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = "rdfhelpers::select({idxs},\n    [{captures}] ( {i} ) {{ return {predExpr}; }})".format(
                    idxs=defCache(self.rng._idxs.op),
                    captures=captures,
                    i="{0} {1}".format(self._i.typeHint, self._i.name),
                    predExpr=defCache(self.predExpr)
                    )
            if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
                return expr
            else:
                funName = defCache.symbol(expr, resultType="ROOT::VecOps::RVec<{0}>".format(SizeType), args=paramDecl)
                return "{0}({1})".format(funName, paramCall)

class Next(TupleOp):
    """ Define a search (first matching item, for a version that processes the whole range see Reduce) """
    def __init__(self, rng, pred):
        self.rng = rng ## PROXY
        self._i = LocalVariablePlaceholder(SizeType)
        self.predExpr = adaptArg(pred(self.rng._base[self._i.result]))
        super(Next, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.predExpr):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp != self._i:
                        yield dp
    @property
    def result(self):
        return self.rng._base[self]
    def __eq__(self, other):
        return isinstance(other, Next) and ( self.rng == other.rng ) and ( self.predExpr == other.predExpr )
    def __repr__(self):
        return "Next({0!r}, {1!r})".format(self.rng, self.predExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.predExpr), (self._i,), defCache=defCache)
        ## should only be a non-empty list for the topmost op in case of nested expressions
        localVarsToName = _collectLocalVars((self.rng, self.predExpr), (self._i,), defCache=defCache)
        with _nameLocalVars(localVarsToName):
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = "rdfhelpers::next({idxs},\n     [{captures}] ( {i} ) {{ return {predexpr}; }})".format(
                    idxs=defCache(self.rng._idxs.op),
                    captures=captures,
                    i="{0} {1}".format(self._i.typeHint, self._i.name),
                    predexpr=defCache(self.predExpr),
                    )
            if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
                return expr
            else:
                funName = defCache.symbol(expr, resultType=SizeType, args=paramDecl)
                return "{0}({1})".format(funName, paramCall)

class Reduce(TupleOp):
    """ Reduce a range to a value (could be a transformation, index...) """
    def __init__(self, rng, start, accuFun):
        self.rng = rng ## PROXY
        self.resultType = start._typeName
        self.start = adaptArg(start)
        self._i = LocalVariablePlaceholder(SizeType)
        self._prevRes = LocalVariablePlaceholder(self.resultType)
        self.accuExpr = adaptArg(accuFun(self._prevRes.result, self.rng._base[self._i.result]))
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in (self.rng, self.start, self.accuExpr):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in (self._i, self._prevRes):
                        yield dp
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.resultType, self)
    def __eq__(self, other):
        return isinstance(other, Reduce) and ( self.rng == other.rng ) and ( self.arg == other.arg ) and ( self.accuExpr == other.accuExpr )
    def __repr__(self):
        return "Reduce({0!r}, {1!r}, {2!r})".format(self.rng, self.start, self.accuExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps((self.rng, self.start, self.accuExpr), (self._i, self._prevRes), defCache=defCache)
        localVarsToName = _collectLocalVars((self.rng, self.start, self.accuExpr), (self._i, self._prevRes), defCache=defCache)
        with _nameLocalVars(localVarsToName):
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = "rdfhelpers::reduce({idxs}, {start},\n     [{captures}] ( {prevRes}, {i} ) {{ return {accuexpr}; }})".format(
                    idxs=defCache(self.rng._idxs.op),
                    start=defCache(self.start),
                    captures=captures,
                    prevRes="{0} {1}".format(self._prevRes.typeHint, self._prevRes.name),
                    i="{0} {1}".format(self._i.typeHint, self._i.name),
                    accuexpr=defCache(self.accuExpr)
                    )
            if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
                return expr
            else:
                funName = defCache.symbol(expr, resultType=self.resultType, args=paramDecl)
                return "{0}({1})".format(funName, paramCall)

class Combine(TupleOp):
    def __init__(self, num, ranges, candPredFun, sameIdxPred=lambda i1,i2: i1 < i2):
        self.n = num
        self.ranges = ranges if len(ranges) > 1 else tuple(repeat(ranges[0], self.n))
        self.candPredFun = candPredFun
        self._i = tuple(LocalVariablePlaceholder(SizeType) for i in range(num))
        from . import treefunctions as op
        areDiff = op.AND(*(sameIdxPred(ia.result, ib.result)
                for ((ia, ra), (ib, rb)) in combinations(zip(self._i, self.ranges), 2)
                if ra._base == rb._base))
        candPred = self.candPredFun(*( rng._base[idx.result] for rng,idx in zip(self.ranges, self._i)))
        if len(areDiff.op.args) > 0:
            self.predExpr = adaptArg(op.AND(areDiff, candPred))
        else:
            self.predExpr = adaptArg(candPred)
    @property
    def resultType(self):
        return "ROOT::VecOps::RVec<rdfhelpers::Combination<{0:d}>>".format(self.n)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if not defCache._getColName(self):
            for arg in chain(self.ranges, [self.predExpr]):
                if select(arg):
                    yield
                for dp in arg.deps(defCache=defCache, select=select, includeLocal=includeLocal):
                    if includeLocal or dp not in self._i:
                        yield dp
    @property
    def result(self):
        from .treeproxies import CombinationListProxy, makeProxy
        return CombinationListProxy(self, makeProxy(self.resultType, self))
    def __eq__(self, other):
        return isinstance(other, Combine) and ( self.n == other.n ) and all( ra == rb for ra,rb in zip(self.ranges, other.ranges) ) and ( self.predExpr == other.predExpr )
    def __repr__(self):
        return "Combine({0:d}, {1!r}, {2!r})".format(self.n, self.ranges, self.predExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList = _collectDeps(chain(self.ranges, [self.predExpr]), self._i, defCache=defCache)
        ## should only be a non-empty list for the topmost op in case of nested expressions
        localVarsToName = _collectLocalVars(chain(self.ranges, [self.predExpr]), self._i, defCache=defCache)
        with _nameLocalVars(localVarsToName):
            captures, paramDecl, paramCall = _convertFunArgs(depList, defCache=defCache)
            expr = ("rdfhelpers::combine{num:d}(\n"
                "     [{captures}] ( {predIdxArgs} ) {{ return {predExpr}; }},\n"
                "     {ranges})").format(
                    num=self.n,
                    captures=captures,
                    predIdxArgs=", ".join("{0} {1}".format(i.typeHint, i.name) for i in self._i),
                    predExpr = defCache(self.predExpr),
                    ranges=", ".join(defCache(rng._idxs.op) for rng in self.ranges)
                    )
            if any(isinstance(dp, LocalVariablePlaceholder) for dp in depList):
                return expr
            else:
                funName = defCache.symbol(expr, resultType=self.resultType, args=paramDecl)
                return "{0}({1})".format(funName, paramCall)

class PsuedoRandom(TupleOp):
    """ Pseudorandom number (integer or float) within range """
    def __init__(self, xMin, xMax, seed, isIntegral=False):
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
