"""
Object representation of operations on TTree branches

The aim is to provide provide sufficiently flexible and complete foundation
for the development of efficient histogram-filling programs
through the use of python wrappers (see e.g. treeproxies).
"""

from itertools import chain

class TupleOp(object):
    """ Interface & base class for operations on leafs and resulting objects / values """
    def deps(self, defCache=None, select=(lambda x : True)):
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        if select(self.length):
            yield self.length
        yield from self.length.deps(defCache=defCache, select=select)
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
    #
    , "abs"   : lambda cppStr,arg : "std::abs( {0} )".format(cppStr(arg))
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
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
        self.index = adaptArg(index, typeHint=SizeType)
        super(GetItem, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (self.arg, self.index):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeName, self)
    def __eq__(self, other):
        return isinstance(other, GetItem) and ( self.arg == other.arg ) and ( self.typeName == other.typeName ) and ( self.index == other.index )
    def __repr__(self):
        return "GetItem({0!r}, {1!r})".format(self.arg, self.index)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        return "{0}[{1}]".format(defCache(self.arg), defCache(self.index))

class Construct(TupleOp):
    def __init__(self, typeName, args):
        self.typeName = typeName
        self.args = tuple(adaptArg(a, typeHint="Double_t") for a in args)
        super(Construct, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
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
        try:
            self._mp  = getattr(ROOT, name)
        except:
            self._mp = None
        self.args = tuple(adaptArg(arg) for arg in args)
        super(CallMethod, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in self.args:
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
    @property
    def result(self):
        retTypeN = next( tok.strip("*&") for tok in self._mp.func_doc.split() if tok != "const" ) if self._mp else "Float_t"
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
        return "{0}({1})".format(self.name, ", ".join(defCache(arg) for arg in self.args))

class CallMemberMethod(TupleOp):
    """ Call a member method """
    def __init__(self, this, name, args):
        self.this = adaptArg(this)
        self.name = name ## NOTE can only be a hardcoded string this way
        self._mp  = getattr(this._typ, name)
        self.args = tuple(adaptArg(arg) for arg in args)
        super(CallMemberMethod, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in chain((self.this,), self.args):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        if select(self.this):
            yield self.this
        yield from self.this.deps(defCache=defCache, select=select)
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
    def __init__(self, typeName, definition):
        self.typeName = typeName
        self.definition = definition
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
        return defCache.symbol(self.definition)

class InitList(TupleOp):
    def __init__(self, typeName, elmType, elms):
        self.typeName = typeName
        self.elms = tuple(adaptArg(e, typeHint=elmType) for e in elms)
        super(InitList, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for elm in self.elms:
            if select(elm):
                yield elm
            yield from elm.deps(defCache=defCache, select=select)
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
    def __init__(self, typeHint, name="<changeme>"):
        self.name = name
        self.typeHint = typeHint
        super(LocalVariablePlaceholder, self).__init__()
    @property
    def result(self):
        from .treeproxies import makeProxy
        return makeProxy(self.typeHint, self)
    def __eq__(self, other): ## TODO ??
        return isinstance(other, LocalVariablePlaceholder) and ( self.name == other.name ) and ( self.typeHint == other.typeHint )
    def __repr__(self):
        return "LocalVariablePlaceholder({0!r}, {1!r})".format(self.name, self.typeHint)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=None):
        return str(self.name)

class Select(TupleOp):
    """ Define a selection on a range """
    def __init__(self, rng, pred):
        self.rng = rng ## PROXY
        self.predExpr = adaptArg(pred(self.rng._base[LocalVariablePlaceholder(SizeType).result]))
        super(Select, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (adaptArg(self.rng), self.predExpr):
            if select(arg):
                yield
            yield from arg.deps(defCache=defCache, select=select)
    @property
    def result(self):
        from .treeproxies import VectorProxy ## FIXME not sure
        return VectorProxy(self, "std::vector<{0}>".format(SizeType))
    def __eq__(self, other):
        return isinstance(other, Select) and ( self.rng == other.rng ) and ( self.predExpr == other.predExpr )
    def __repr__(self):
        return "Select({0!r}, {1!r})".format(self.rng, self.predExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList = set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn) or ( defCache.backend.shouldDefine(op) and defCache._getName(op) ))) ## FIXME define as "isLeaf" or so
            for expr in (self.rng, self.predExpr)))
        funName = defCache.symbol((
            "using namespace ROOT::VecOps;\n"
            "RVec<std::size_t> <<name>>({fargs})\n"
            "{{\n"
            "  return rdfhelpers::select({idxs},\n"
            "     [{captures}] ( std::size_t i ) {{ return {predexpr}; }});\n"
            "}};\n"
            ).format(
                fargs=", ".join( ## TODO this could be improved/factored out (grouped with captures)
                    ("const RVec<{0}>& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetArrayColumn)
                    else "const {0}& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetColumn)
                    else "(problem with {0!r})".format(ld)) for ld in depList),
                idxs=defCache(self.rng._idxs.op),
                captures=",".join("&{0}".format(ld.name) for ld in depList),
                predexpr=defCache(self.predExpr).replace("<changeme>", "i")
            ))
        return "{0}({1})".format(funName, ", ".join(ld.name for ld in depList))

class Next(TupleOp):
    """ Define a search (first matching item, for a version that processes the whole range see Reduce) """
    def __init__(self, rng, pred):
        self.rng = rng ## PROXY
        self._i = LocalVariablePlaceholder(SizeType).result
        self.predExpr = adaptArg(pred(self.rng._base[self._i]))
        super(Next, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (self.rng, self.predExpr):
            if select(arg):
                yield
            yield from arg.deps(defCache=defCache, select=select)
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
        depList = set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn) or ( defCache.backend.shouldDefine(op) and defCache._getName(op) ))) ## FIXME define as "isLeaf" or so
            for expr in (self.rng, self.predExpr)))
        funName = defCache.symbol((
            "using namespace ROOT::VecOps;\n"
            "std::size_t <<name>>({fargs})\n"
            "{{\n"
            "  return rdfhelpers::next({idxs},\n"
            "     [{captures}] ( std::size_t i ) {{ return {predexpr}; }});\n"
            "}};\n"
            ).format(
                fargs=", ".join( ## TODO this could be improved/factored out (grouped with captures)
                    ("const RVec<{0}>& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetArrayColumn)
                    else "const {0}& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetColumn)
                    else "(problem with {0!r})".format(ld)) for ld in depList),
                idxs=defCache(self.rng._idxs.op),
                captures=",".join("&{0}".format(ld.name) for ld in depList),
                predexpr=defCache(self.predExpr).replace("<changeme>", "i")
            ))
        return "{0}({1})".format(funName, ", ".join(ld.name for ld in depList))

class Reduce(TupleOp):
    """ Reduce a range to a value (could be a transformation, index...) """
    def __init__(self, rng, start, accuFun):
        self.rng = rng ## PROXY
        self.resultType = start._typeName
        self.start = adaptArg(start)
        self._i = LocalVariablePlaceholder(SizeType, name="<changeme_index>").result
        self._prevRes = LocalVariablePlaceholder(self.resultType, name="<changeme_result>").result
        self.accuExpr = adaptArg(accuFun(self._prevRes, self.rng._base[self._i]))
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (self.rng, self.start, self.accuExpr):
            if select(arg):
                yield
            yield from arg.deps(defCache=defCache, select=select)
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
        depList = set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn))) ## FIXME define as "isLeaf" or so
            for expr in (self.rng, self.start, self.accuExpr)))
        from .treeproxies import ListBase
        for expr in (self.rng, self.start, self.accuExpr):
            for dep in expr.deps(defCache=defCache, select=(lambda op : defCache.backend.shouldDefine(op) and defCache._getColName(op))):
                depResult = dep.result
                if isinstance(depResult, ListBase):
                    depList.add(GetArrayColumn(depResult.valueType, defCache._getColName(dep), adaptArg(depResult.__len__())))
                else:
                    depList.add(GetColumn(depResult._typeName, defCache._getColName(dep)))
        funName = defCache.symbol((
            "using namespace ROOT::VecOps;\n"
            "{resType} <<name>>({fargs})\n"
            "{{\n"
            "  return rdfhelpers::reduce({idxs}, {start},\n"
            "     [{captures}] ( {resType} prevResult, std::size_t i ) {{ return {accuexpr}; }});\n"
            "}};\n"
            ).format(
                resType=self.resultType,
                fargs=", ".join( ## TODO this could be improved/factored out (grouped with captures)
                    ("const RVec<{0}>& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetArrayColumn)
                    else "const {0}& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetColumn)
                    else "(problem with {0!r})".format(ld)) for ld in depList),
                idxs=defCache(self.rng._idxs.op),
                start=defCache(self.start),
                captures=",".join("&{0}".format(ld.name) for ld in depList),
                accuexpr=defCache(self.accuExpr).replace("<changeme_index>", "i").replace("<changeme_result>", "prevResult")
            ))
        return "{0}({1})".format(funName, ", ".join(ld.name for ld in depList))

class KinematicVariation(TupleOp):
    def __init__(self, rng, modif, pred):
        self.rng = rng ## PROXY (original range)
        self.modifExpr = adaptArg(modif(self.rng._base[LocalVariablePlaceholder(SizeType).result]))
        p4local = LocalVariablePlaceholder("ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>")
        p4local.name = "<modp4>"
        self.predExpr = adaptArg(pred(p4local.result, self.rng._base[LocalVariablePlaceholder(SizeType).result]))
        super(KinematicVariation, self).__init__()
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (adaptArg(self.rng), self.modifExpr, self.predExpr):
            if select(arg):
                yield
            yield from arg.deps(defCache=defCache, select=select)
    @property
    def result(self):
        from .treeproxies import makeProxy, ModifiedCollectionProxy
        return ModifiedCollectionProxy(makeProxy("rdfhelpers::ModifiedKinCollection", self), self.rng._base) ## itemType to be filled in by 
    def __eq__(self, other):
        return isinstance(other, KinematicVariation) and ( self.rng == other.rng ) and ( self.modifExpr == other.modifExpr ) and ( self.predExpr == other.predExpr )
    def __repr__(self):
        return "KinematicVariation({0!r}, {1!r}, {2!r})".format(self.rng, self.modifExpr, self.predExpr)
    def __hash__(self):
        return hash(self.__repr__())
    def get_cppStr(self, defCache=cppNoRedir):
        depList_modif = set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn) or ( defCache.backend.shouldDefine(op) and defCache._getName(op) ))) ## FIXME define as "isLeaf" or so
            for expr in (self.rng, self.modifExpr)))
        depList_pred = set(chain.from_iterable(
            expr.deps(defCache=defCache, select=(lambda op : isinstance(op, GetColumn) or isinstance(op, GetArrayColumn) or ( defCache.backend.shouldDefine(op) and defCache._getName(op) ))) ## FIXME define as "isLeaf" or so ==> can have a 'Lambda' op taking care of this
            for expr in (self.rng, self.predExpr)))
        depList_merged = set(depList_modif)
        depList_merged.update(depList_pred)
        funName = defCache.symbol((
            "using namespace ROOT::VecOps;\n"
            "rdfhelpers::ModifiedKinCollection <<name>>({fargs})\n"
            "{{\n"
            "  return rdfhelpers::modifyKinCollection({idxs},\n"
            "     [{captures_modif}] ( std::size_t i ) {{ return {modifExpr}; }},\n"
            "     [{captures_pred}] ( const rdfhelpers::ModifiedKinCollection::LorentzVector& p4Mod, std::size_t i ) {{ return {predExpr}; }}\n"
            "  );\n"
            "}};\n"
            ).format(
                fargs=", ".join( ## TODO this could be improved/factored out (grouped with captures)
                    ("const RVec<{0}>& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetArrayColumn)
                    else "const {0}& {1}".format(ld.typeName, ld.name) if isinstance(ld, GetColumn)
                    else "(problem with {0!r})".format(ld)) for ld in depList_merged),
                idxs=defCache(self.rng._idxs.op),
                captures_modif=",".join("&{0}".format(ld.name) for ld in depList_modif),
                captures_pred=",".join("&{0}".format(ld.name) for ld in depList_pred),
                modifExpr=defCache(self.modifExpr).replace("<changeme>", "i"),
                predExpr=defCache(self.predExpr).replace("<changeme>", "i").replace("<modp4>", "p4Mod")
            ))
        return "{0}({1})".format(funName, ", ".join(ld.name for ld in depList_merged))
