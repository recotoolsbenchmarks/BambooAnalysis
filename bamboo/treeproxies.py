"""
Wrapper classes for TTree/ROOT::RDataFrame

The goal is to make it easier to work with flat trees by generating structs
"on the fly", and allow postponing evaluation to a later stage
(as in RDataFrame, or previously with a code-generating backend).
This is done by providing proxy objects to parts of the record content,
and "operations" derived from these, which are again wrapped in proxies,
to allow easy constructing the expression tree for the backend.
As an example: 'myTreeProxy.event_id' will return an integer proxy,
with as 'op' attribute a 'GetColumn(tree, "event_id")' operation.
"""

from .treeoperations import *
import functools

boolType = "bool"
_boolTypes = set(("bool", "Bool_t"))
intType = "Int_t"
_integerTypes = set(("Int_t", "UInt_t", "Char_t", "UChar_t", "ULong64_t", "int", "unsigned", "unsigned short", "char", "signed char", "unsigned char", "Short_t", "size_t", "std::size_t", "unsigned short", "unsigned long"))
floatType = "Double_t"
_floatTypes = set(("Float_t", "Double_t", "float", "double"))
## TODO lists are probably not complete
import re
vecPat = re.compile(r"(?:vector|ROOT\:\:VecOps\:\:RVec)\<(?P<item>[a-zA-Z_0-9\<\>,\: ]+)\>")

def makeProxy(typeName, parent, length=None):
    if length is not None:
        return ArrayProxy(parent, typeName, length)
    if typeName in _boolTypes:
        return BoolProxy(parent, typeName)
    elif typeName in _integerTypes:
        return IntProxy(parent, typeName)
    elif typeName in _floatTypes:
        return FloatProxy(parent, typeName)
    else:
        m = vecPat.match(typeName)
        if m:
            return VectorProxy(parent, typeName)
        else:
            return ObjectProxy(parent, typeName)
def makeConst(value, typeHint):
    return makeProxy(typeHint, adaptArg(value, typeHint))

class NumberProxy(TupleBaseProxy):
    def __init__(self, parent, typeName):
        super(NumberProxy, self).__init__(typeName, parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._typeName)
    def _binaryOp(self, opName, other, outType=None):
        if outType is None:
            outType = self._typeName ## default: math will give same type
        return MathOp(opName, self, other, outType=outType).result
for nm,opNm in {
          "__add__": "add"
        , "__sub__": "subtract"
        , "__mul__": "multiply"
        , "__pow__": "pow"
        }.items():
    setattr(NumberProxy, nm, functools.partialmethod(
        (lambda self, oN, other : self._binaryOp(oN, other)), opNm))
for nm in ("__lt__", "__le__", "__eq__", "__ne__", "__gt__", "__ge__"):
    setattr(NumberProxy, nm, (lambda oN, oT : (lambda self, other : self._binaryOp(oN, other, outType=oT)))(nm.strip("_"), boolType))

class IntProxy(NumberProxy):
    def __init__(self, parent, typeName):
        super(IntProxy, self).__init__(parent, typeName)
for nm,(opNm,outType) in {
          "__truediv__" : ("floatdiv", floatType)
        , "__floordiv__": ("divide"  , None)
        , "__and__"     : ("band"    , None)
        , "__or__"      : ("bor"     , None)
        , "__invert__"  : ("bxor"    , None)
        , "__xor__"     : ("bnot"    , None)
        }.items():
    setattr(IntProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, outType))

class BoolProxy(IntProxy):
    """ Proxy for a boolean type """
    def __init__(self, parent, typeName):
        super(BoolProxy, self).__init__(parent, typeName)
    def __repr__(self):
        return "BoolProxy({0!r}, {1!r})".format(self._parent, self._typeName)
for nm,opNm in {
          "__and__"   : "and"
        , "__or__"    : "or"
        , "__invert__": "not"
        , "__xor__"   : "ne"
        }.items():
    setattr(BoolProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, boolType))

class FloatProxy(NumberProxy):
    def __init__(self, parent, typeName):
        super(FloatProxy, self).__init__(parent, typeName)
for nm,(opNm,outType) in {
          "__truediv__": ("div" , None)
        }.items():
    setattr(FloatProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, outType))

def _hasFloat(*args):
    return any(isinstance(a, float) or isinstance(a, FloatProxy) for a in args)

class ArrayProxy(TupleBaseProxy):
    """ (possibly var-sized) array of anything """
    def __init__(self, parent, typeName, length):
        self._length = length
        super(ArrayProxy, self).__init__("{0}[]".format(typeName), parent=parent)
        self.valueType = typeName
    def __getitem__(self, index):
        return GetItem(self, self.valueType, index).result
    def __len__(self):
        return self._length.result
    def __repr__(self):
        return "ArrayProxy({0!r}, {1!r}, {2!r})".format(self._parent, self._typeName, self._length)

class LeafGroupProxy(TupleBaseProxy):
    """ Base class for proxies with a prefix (leaf group, container) """
    def __init__(self, prefix, parent):
        self._prefix = prefix
        super(LeafGroupProxy, self).__init__("struct", parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._prefix, self._parent)
def fullPrefix(prefix, parent): ## TODO check if we actually need it (probably yes for nested)
    if parent:
        return "".join(fullPrefix(parent._prefix, parent._parent), prefix)
    else:
        return prefix

class TreeBaseProxy(LeafGroupProxy):
    """ Tree proxy base class """
    def __init__(self, tree):
        self._tree = tree
        super(TreeBaseProxy, self).__init__("", None)
    @property
    def op(self):
        return None
    def __repr__(self):
        return "{0}({1!r})".format(self.__class__.__name__, self._tree)
    def deps(self, defCache=None, select=(lambda x : True), includeLocal=False):
        yield from []

class ListBase(object):
    """ Interface definition for range proxies (Array/Vector, split object vector, selection/reordering)

    Important for users: _base always contains a basic container (e.g. ContainerGroupProxy), and _idxs
    the indices out of _base this list contains (so _base[_idxs[i]] for i in range(len) are always valid).

    TODO verify / stress-tests (selection of selection, next of selection of selection etc.)
    """
    def __init__(self):
        self.valueType = None
        self._base = self ## TODO get rid of _
        super(ListBase, self).__init__()
    def __getitem__(self, index):
        pass ## need override
    def __len__(self):
        pass ## need overridde
    @property
    def _idxs(self):
        return Construct("rdfhelpers::IndexRange<{0}>".format(SizeType), (adaptArg(self.__len__()),)).result ## FIXME uint->int narrowing
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from []

class ContainerGroupItemProxy(TupleBaseProxy):
    """ Proxy for an item in a structure of arrays """
    def __init__(self, parent, idx):
        self._idx = adaptArg(idx, typeHint=SizeType)
        super(ContainerGroupItemProxy, self).__init__("struct", parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._idx)

class ContainerGroupProxy(LeafGroupProxy,ListBase):
    """ Proxy for a structure of arrays """
    def __init__(self, prefix, parent, size, valuetype):
        ListBase.__init__(self)
        self._size = size
        self.valuetype = valuetype
        super(ContainerGroupProxy, self).__init__(prefix, parent)
    def __len__(self):
        return self._size.result
    def __getitem__(self, index):
        return self.valuetype(self, index)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._size)
    ##
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from self._size.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self, other):
        return isinstance(other, ContainerGroupProxy) and ( self._size == other._size ) and ( self.valuetype == other.valuetype ) and ( self._parent == other._parent ) and ( self._prefix == other._prefix )

class ObjectProxy(NumberProxy):
    """
    Imitate an object
    """
    __slots__ = ("_typ",)
    def __init__(self, parent, typeName):
        from cppyy import gbl
        self._typ = getattr(gbl, typeName) ## NOTE could also use TClass machinery
        super(ObjectProxy, self).__init__(parent, typeName)
    def __getattr__(self, name):
        if name not in dir(self._typ):
            raise AttributeError("Type {0} has no member {1}".format(self._typeName, name))
        from cppyy import gbl
        if hasattr(self._typ, name) and isinstance(getattr(self._typ, name), gbl.MethodProxy):
            return ObjectMethodProxy(self, name)
        else:
            return GetDataMember(self, name).result
    def __repr__(self):
        return "ObjectProxy({0!r}, {1!r})".format(self._parent, self._typeName)

class ObjectMethodProxy(TupleBaseProxy): ## TODO data members?
    """
    Imitate a member method of an object
    """
    __slots__ = ("_objStb", "_name")
    def __init__(self, objStb, name):
        self._objStb = objStb
        self._name = name
        super(ObjectMethodProxy, self).__init__("{0}.{1}(...)".format(objStb._typeName, self._name))
    def __call__(self, *args):
        ## TODO maybe this is a good place to resolve the right overload? or do some arguments checking
        return CallMemberMethod(self._objStb, self._name, tuple(args)).result
    def __repr__(self):
        return "ObjectMethodProxy({0!r}, {1!r})".format(self._objStb, self._name)

class MethodProxy(TupleBaseProxy):
    """
    Imitate a free-standing method
    """
    __slots__ = ("_name",)
    def __init__(self, name):
        self._name = name
        super(MethodProxy, self).__init__("{0}(...)".format(self._name))
    def __call__(self, *args):
        ## TODO maybe this is a good place to resolve the right overload? or do some arguments checking
        return CallMethod(self._name, tuple(args)).result
    def __repr__(self):
        return "MethodProxy({0!r})".format(self._name)

class VectorProxy(ObjectProxy,ListBase):
    """ Vector-as-array (maybe to be eliminated with var-array branches / generalised into object) """
    def __init__(self, parent, typeName):
        ListBase.__init__(self)
        from cppyy import gbl
        vecClass = getattr(gbl, typeName)
        if hasattr(vecClass, "value_type"):
            value = getattr(vecClass, "value_type")
            if hasattr(value, "__cppname__"):
                self.valueType = value.__cppname__
            elif str(value) == value:
                self.valueType = value
            else:
                raise RuntimeError("value_type attribute of {0} is neither a PyROOT class nor a string, but {1}".format(typeName, value))
        else:
            self.valueType = vecPat.match(typeName).group("item")
        super(VectorProxy, self).__init__(parent, typeName)
    def __getitem__(self, index):
        return GetItem(self, self.valueType, index).result
    def __len__(self):
        return self.size()
    def __repr__(self):
        return "VectorProxy({0!r}, {1!r})".format(self._parent, self._typeName)
    #
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from self.parent.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self, other):
        return isinstance(other, VectorProxy) and ( self._parent == other._parent ) and ( self.valuetype == other.valuetype )

class SelectionProxy(TupleBaseProxy,ListBase):
    """ Proxy for a selection from an iterable (ContainerGroup/ other selection etc.) """
    def __init__(self, parent):
        ListBase.__init__(self)
        self._base = parent.rng._base
        self.valueType = self._base.valueType
        ## the list of indices is stored as the parent
        super(SelectionProxy, self).__init__(self.valueType, parent=parent)
    @property
    def _idxs(self):
        return self._parent.result
    def __getitem__(self, index):
        return self._base[GetItem(self._idxs, SizeType, index)]
    def __len__(self):
        return self._idxs.__len__()
    def __repr__(self):
        return "SelectionProxy({0!r})".format(self._parent)
    #
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if select(self._parent):
            yield self._parent
        yield from self._parent.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self, other):
        return isinstance(other, SelectionProxy) and ( self._parent == other._parent ) and ( self.valueType == other.valueType ) and ( self._base == other._base )

## Attribute classes (like property) to customize the above base classes

class proxy(object): ## the default one
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op.result

class funProxy(object): ## the generic one
    def __init__(self, fun):
        self.fun = fun
    def __get__(self, inst, cls):
        return self.fun(inst)

class itemProxy(object):
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op[inst._idx]

class itemRefProxy(object):
    def __init__(self, op, getTarget):
        self.op = op
        self.getTarget = getTarget
    def __get__(self, inst, cls):
        return self.getTarget(inst)[self.op[inst._idx]]

class itemObjProxy(object): ## re-construct an object that was split in arrays
    def __init__(self, typeName, args):
        self.typeName = typeName
        self.args = args
    def __get__(self, inst, cls):
        return Construct(self.typeName, tuple(arg[inst._idx] for arg in self.args)).result

class varItemProxy(object):
    def __init__(self, op):
        self.op = op
    def __get__(self, inst, cls):
        return self.op[inst._parent._parent.result.indices()[inst._idx]]
class varItemRefProxy(object):
    def __init__(self, op, getTarget):
        self.op = op
        self.getTarget = getTarget
    def __get__(self, inst, cls):
        return self.getTarget(inst)[self.op[inst._parent._parent.result.indices()[inst._idx]]]

class JMEVariations(TupleBaseProxy):
    def __init__(self, parent, orig, args, varItemType=None):
        self.orig = orig
        self._args = args
        self.varItemType = varItemType if varItemType else orig.valuetype
        self.calc = None
        self.calcProd = DefinedVar("JMESystematicsCalculator", "JMESystematicsCalculator <<name>>;").result.produceModifiedCollections(*self._args)
    def initCalc(self, calcProxy, calcHandle=None):
        self.calc = calcHandle ## handle to the actual module object
        self.calcProd = calcProxy.produceModifiedCollections(*self._args)
    def __getitem__(self, key):
        if self.calc and not self.calc.hasProduct(key):
            raise KeyError("Modified jet collection with name {0!r} will not be produced, please check the configuration".format(key))
        #res_item = GetItem(self.calcProd, "JMESystematicsCalculator::result_entry_t", makeConst(key, "std::string"), "std::string")
        res_item = self.calcProd.at(makeConst(key, "std::string")).op
        return ModifiedCollectionProxy(res_item, self.orig, itemType=self.varItemType)

class ModifiedCollectionProxy(TupleBaseProxy,ListBase):
    ## TODO make jet-specific, maybe add a bit of MET deco
    """ Collection with a selection and modified branches """
    def __init__(self, parent, base, itemType=None):
        ## parent has a member indices() and one for each of modBranches
        ListBase.__init__(self)
        self.orig = base
        self.valueType = base.valueType ## for ListBase
        self.itemType = itemType ## for actual items
        super(ModifiedCollectionProxy, self).__init__(self.valueType, parent=parent)
    def __len__(self):
        return self._parent.result.indices().__len__()
    def __getitem__(self, index):
        return self.itemType(self, index)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self.itemType)
    ##
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        if select(self._parent):
            yield self._parent
        yield from self._parent.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self, other):
        return isinstance(other, self.__class__) and ( self._parent == other._parent ) and ( self.itemType == other.itemType )

class CombinationProxy(TupleBaseProxy):
    ## NOTE decorated rdfhelpers::Combination
    def __init__(self, cont, parent):
        self.cont = cont
        super(CombinationProxy, self).__init__("struct", parent=parent) ## parent=ObjectProxy for a combination (indices)
    @property
    def _idx(self):
        return adaptArg(self._parent)._index ## parent is a GetIndex-proxy
    @property
    def index(self):
        return adaptArg(self._parent).index ## parent is a GetIndex-proxy
    def __getitem__(self, i):
        idx = makeConst(i, SizeType)
        return self.cont.base(i)[self._parent.get(idx)]
    ## TODO add more (maybe simply defer to rdfhelpers::Combination object
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.cont, self._parent)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from self._parent.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self):
        return isinstance(other, self.__class__) and self._parent == other._parent

class CombinationListProxy(TupleBaseProxy,ListBase):
    ## NOTE list of decorated rdfhelpers::Combination
    ## TODO check if this works with selection (it should...)
    def __init__(self, parent, combs):
        self.combs = combs ## straightforward proxy for parent
        ListBase.__init__(self) ## FIXME check above how to do this correctly...
        super(CombinationListProxy, self).__init__(parent.resultType, parent=parent)
    def __getitem__(self, idx):
        return CombinationProxy(self, GetItem(self.combs, self.combs.valueType, adaptArg(idx, SizeType)).result)
    def base(self, i):
        return self._parent.ranges[i]._base
    def __len__(self):
        return self.combs.size()
    def __repr__(self):
        return "{0}({1!r})".format(self.__class__.__name__, self._parent)
    def deps(self, defCache=cppNoRedir, select=(lambda x : True), includeLocal=False):
        yield from self._parent.deps(defCache=defCache, select=select, includeLocal=includeLocal)
    def __eq__(self):
        return isinstance(other, self.__class__) and self._parent == other._parent
