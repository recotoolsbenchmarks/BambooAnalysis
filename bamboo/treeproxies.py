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
            return VectorProxy(parent, typeName=typeName)
        else:
            return ObjectProxy(parent, typeName)
def makeConst(value, typeHint):
    return makeProxy(typeHint, adaptArg(value, typeHint))

class NumberProxy(TupleBaseProxy):
    def __init__(self, parent, typeName):
        super(NumberProxy, self).__init__(typeName, parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._typeName)
    def _binaryOp(self, opName, other, outType=None, reverseOrder=False):
        if outType is None:
            outType = self._typeName ## default: math will give same type
        if reverseOrder:
            return MathOp(opName, other, self, outType=outType).result
        return MathOp(opName, self, other, outType=outType).result
    def _unaryOp(self, opName):
        return MathOp(opName, self, outType=self._typeName)
for nm,opNm in {
          "__add__": "add"
        , "__sub__": "subtract"
        , "__mul__": "multiply"
        , "__pow__": "pow"
        , "__radd__": "add"
        , "__rmul__": "multiply"
        }.items():
    setattr(NumberProxy, nm, functools.partialmethod(
        (lambda self, oN, other : self._binaryOp(oN, other)), opNm))
for nm,opNm in {
          "__rsub__": "subtract"
        , "__rpow__": "rpow"
        }.items():
    setattr(NumberProxy, nm, functools.partialmethod(
        (lambda self, oN, other : self._binaryOp(oN, other, reverseOrder=True)), opNm))
for nm in ("__lt__", "__le__", "__eq__", "__ne__", "__gt__", "__ge__"):
    setattr(NumberProxy, nm, (lambda oN, oT : (lambda self, other : self._binaryOp(oN, other, outType=oT)))(nm.strip("_"), boolType))
for name in ("__neg__",):
    setattr(NumberProxy, name, functools.partialmethod(
        (lambda self, oN: self._unaryOp(oN)), name.strip("_")))

class IntProxy(NumberProxy):
    def __init__(self, parent, typeName):
        super(IntProxy, self).__init__(parent, typeName)
for nm,(opNm,outType) in {
          "__truediv__" : ("floatdiv", floatType)
        , "__floordiv__": ("divide"  , None)
        , "__mod__"     : ("mod"     , None)
        , "__lshift__"  : ("lshift"  , None)
        , "__rshift__"  : ("rshift"  , None)
        , "__and__"     : ("band"    , None)
        , "__or__"      : ("bor"     , None)
        , "__xor__"     : ("bxor"    , None)
        , "__rand__"     : ("band"    , None)
        , "__ror__"      : ("bor"     , None)
        , "__rxor__"     : ("bxor"    , None)
        }.items():
    setattr(IntProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, outType))
for nm,(opNm,outType) in {
          "__rtruediv__" : ("floatdiv", floatType)
        , "__rfloordiv__": ("divide"  , None)
        , "__rmod__"     : ("mod"     , None)
        , "__rlshift__"  : ("lshift"  , None)
        , "__rrshift__"  : ("rshift"  , None)
        }.items():
    setattr(IntProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT, reverseOrder=True)),
        opNm, outType))
for name,opName in {
            "__invert__": "bnot"
        }.items():
    setattr(IntProxy, name, functools.partialmethod(
        (lambda self, oN: self._unaryOp(oN)), opName))

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
        , "__rand__"   : "and"
        , "__ror__"    : "or"
        , "__rxor__"   : "ne"
        }.items():
    setattr(BoolProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, boolType))
for name,opName in {
            "__invert__": "not"
        }.items():
    setattr(BoolProxy, name, functools.partialmethod(
        (lambda self, oN: self._unaryOp(oN)), opName))

class FloatProxy(NumberProxy):
    def __init__(self, parent, typeName):
        super(FloatProxy, self).__init__(parent, typeName)
for nm,(opNm,outType) in {
          "__truediv__": ("floatdiv" , None)
        }.items():
    setattr(FloatProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT)),
        opNm, outType))
for nm,(opNm,outType) in {
          "__rtruediv__": ("floatdiv" , None)
        }.items():
    setattr(FloatProxy, nm, functools.partialmethod(
        (lambda self, oN, oT, other : self._binaryOp(oN, other, outType=oT, reverseOrder=True)),
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
        if isinstance(index, slice):
            if index.step and index.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self, index.start, index.stop, valueType=self.valueType)
        else:
            return GetItem(self, self.valueType, index).result
    def __len__(self):
        return self._length
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
    def __eq__(self, other):
        return self._tree == other._tree
    def __deepcopy__(self, memo):
        ## *never* copy the TTree, although copying proxies is fine
        return self.__class__(self._tree)

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
        """ Get item using index of this range """
        if isinstance(index, slice):
            if index.step and index.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self, index.start, index.stop, valueType=self.valueType)
        else:
            return self._getItem(index)
    def _getItem(self, baseIndex):
        """ Get item using index of base range """
        return self._base[baseIndex]
    def __len__(self):
        pass ## need overridde
    @property
    def _idxs(self):
        return Construct("rdfhelpers::IndexRange<{0}>".format(SizeType), (adaptArg(self.__len__()),)).result ## FIXME uint->int narrowing

class ContainerGroupItemProxy(TupleBaseProxy):
    """ Proxy for an item in a structure of arrays """
    def __init__(self, parent, idx):
        self._idx = adaptArg(idx, typeHint=SizeType)
        super(ContainerGroupItemProxy, self).__init__("struct", parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._idx)
    @property
    def isValid(self):
        return self._idx.result != -1
    @property
    def idx(self):
        return self._idx.result
    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            raise RuntimeError("Cannot compare proxies of different types: {0} and {1}".format(other.__class__.__name__, self.__class__.__name__))
        if self._parent._base != other._parent._base:
            raise RuntimeError("Cannot compare elements from different containers")
        return self._idx.result == other._idx.result
    def __ne__(self, other):
        from bamboo.treefunctions import NOT
        return NOT(self == other)

class ContainerGroupProxy(LeafGroupProxy,ListBase):
    """ Proxy for a structure of arrays """
    def __init__(self, prefix, parent, size, valuetype):
        ListBase.__init__(self)
        self._size = size
        self.valuetype = valuetype
        super(ContainerGroupProxy, self).__init__(prefix, parent)
    def __len__(self):
        return self._size.result
    def _getItem(self, baseIndex):
        return self.valuetype(self, baseIndex)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self._size)

class ObjectProxy(NumberProxy):
    """
    Imitate an object
    """
    __slots__ = tuple()
    def __init__(self, parent, typeName):
        super(ObjectProxy, self).__init__(parent, typeName)
    @property
    def _typ(self):
        from .root import gbl
        return getattr(gbl, self._typeName)
    def __getattr__(self, name):
        typ = self._typ
        if name not in dir(typ):
            raise AttributeError("Type {0} has no member {1}".format(self._typeName, name))
        from .root import gbl
        if hasattr(typ, name) and isinstance(getattr(typ, name), gbl.MethodProxy):
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
    __slots__ = ("_name", "_returnType", "_getFromRoot")
    def __init__(self, name, returnType=None, getFromRoot=True):
        self._name = name
        self._returnType = returnType
        self._getFromRoot = getFromRoot
        super(MethodProxy, self).__init__("{0}(...)".format(self._name))
    def __call__(self, *args):
        ## TODO maybe this is a good place to resolve the right overload? or do some arguments checking
        return CallMethod(self._name, tuple(args), returnType=self._returnType, getFromRoot=self._getFromRoot).result
    def __repr__(self):
        return "MethodProxy({0!r})".format(self._name)

class VectorProxy(ObjectProxy,ListBase):
    """ Vector-as-array (maybe to be eliminated with var-array branches / generalised into object) """
    def __init__(self, parent, typeName=None, itemType=None):
        ListBase.__init__(self)
        if itemType:
            self.valueType = itemType
        else:
            from .root import gbl
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
    def _getItem(self, index):
        return GetItem(self, self.valueType, index).result
    def __len__(self):
        return self.size()
    def __repr__(self):
        return "VectorProxy({0!r}, {1!r})".format(self._parent, self._typeName)

class SelectionProxy(TupleBaseProxy,ListBase):
    """ Proxy for a selection from an iterable (ContainerGroup/ other selection etc.) """
    def __init__(self, base, parent, valueType=None): ## 'parent' is a Select or Sort TupleOp
        ListBase.__init__(self)
        self._base = base
        self.valueType = valueType if valueType else self._base.valueType
        ## the list of indices is stored as the parent
        super(SelectionProxy, self).__init__(self.valueType, parent=parent)
    @property
    def _idxs(self):
        return self._parent.result
    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step and index.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self, index.start, index.stop, valueType=self.valueType)
        else:
            return self._getItem(self._idxs[index])
    def _getItem(self, baseIndex):
        itm = self._base[baseIndex]
        if self.valueType and self.valueType != self._base.valueType:
            return self.valueType(itm._parent, itm._idx)
        else:
            return itm
    def __len__(self):
        return self._idxs.__len__()
    def __repr__(self):
        return "SelectionProxy({0!r}, {1!r}, valueType={2!r})".format(self._base, self._parent, self.valueType)

class SliceProxy(TupleBaseProxy,ListBase):
    """ Proxy for part of an iterable (ContainerGroup/selection etc.) """
    def __init__(self, parent, begin, end, valueType=None): ## 'parent' is another proxy (ListBase, will become _base of this one)
        ListBase.__init__(self)
        self._base = parent if parent is not None else None
        self._begin = makeProxy(SizeType, adaptArg(begin, SizeType)) if begin is not None else None ## None signals 0
        self._end = makeProxy(SizeType, adaptArg(end, SizeType)) if end is not None else makeConst(parent.__len__(), SizeType)
        self.valueType = valueType if valueType else parent.valueType
        super(SliceProxy, self).__init__(self.valueType, parent=adaptArg(parent))
    def _offset(self, idx):
        if self._begin is not None:
            return self._begin+idx
        else:
            return idx
    @property
    def begin(self):
        if self._begin:
            return self._begin
        else:
            return makeConst(0, SizeType)
    @property
    def _idxs(self):
        return Construct("rdfhelpers::IndexRange<{0}>".format(SizeType), (
            adaptArg(self.begin), adaptArg(self._end))).result ## FIXME uint->int narrowing
    def __getitem__(self, key):
        if isinstance(key, slice):
            if key.step and key.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self._base,
                    (self._offset(key.start) if key.start is not None else self._begin),
                    (self._offset(key.stop ) if key.stop  is not None else self._end  ),
                    valueType=self.valueType)
        else:
            return self._getItem(self._offset(key))
    def _getItem(self, baseIndex):
        itm = self._base[baseIndex]
        if self.valueType and self.valueType != self._base.valueType:
            return self.valueType(itm._parent, itm._idx)
        return itm
    def __len__(self):
        return self._end-self.begin
    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r}, valueType={4!r})".format(self.__class__.__name__, self._parent, self._begin, self._end, self.valueType)

class AltLeafVariations(TupleBaseProxy):
    """ Branch with alternative names """
    def __init__(self, parent, brMap, typeName):
        super(AltLeafVariations, self).__init__(typeName, parent=parent)
        self.brMap = brMap
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise ValueError("Getting variation with non-string key {0!r}".format(key))
        if key not in self.brMap:
            raise KeyError("No such variation: {0}".format(key))
        return self.brMap[key].result

class AltLeafGroupVariations(TupleBaseProxy):
    """ Set of groups with alternative names for some branches """
    def __init__(self, parent, orig, brMapMap, altProxyType):
        self.orig = orig
        self.brMapMap = brMapMap
        self.altProxyType = altProxyType
        super(AltLeafGroupVariations, self).__init__("AltLeafGroupVariations", parent=parent)
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise ValueError("Getting variation with non-string key {0!r}".format(key))
        if key not in self.brMapMap:
            raise KeyError("No such variation: {0}".format(key))
        return self.altProxyType(self._parent, self.orig, self.brMapMap[key])

class AltLeafGroupProxy(TupleBaseProxy):
    """ Group with alternative names for some branches """
    ## base class like LeafGroupProxy, but with a brMap
    def __init__(self, parent, orig, brMap, typeName="struct"):
        self.orig = orig
        self.brMap = brMap
        super(AltLeafGroupProxy, self).__init__(typeName, parent=parent)
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self._parent, self.brMap)

class AltCollectionVariations(TupleBaseProxy):
    """ Set of collections with alternative names for some branches """
    def __init__(self, parent, orig, brMapMap, altItemType=None):
        self.orig = orig
        self.brMapMap = brMapMap
        self.altItemType = altItemType if altItemType else orig.valuetype
        super(AltCollectionVariations, self).__init__("AltCollectionVariations", parent=parent)
    def __getitem__(self, key):
        if not isinstance(key, str):
            raise ValueError("Getting variation with non-string key {0!r}".format(key))
        if key not in self.brMapMap:
            raise KeyError("No such variation: {0}. If stored on the NanoAOD: are the branch names in the correct format; if calculated on the fly: has the calculator been configured?".format(key))
        return AltCollectionProxy(self._parent, self.orig, self.brMapMap[key], itemType=self.altItemType)

class AltCollectionProxy(TupleBaseProxy, ListBase):
    """ Collection with alternative names for some branches """
    def __init__(self, parent, orig, brMap, itemType=None):
        ListBase.__init__(self)
        self.orig = orig
        self.valueType = orig.valueType ## for ListBase
        self.itemType = itemType ## for actual items
        self.brMap = brMap
        super(AltCollectionProxy, self).__init__(self.valueType, parent=parent)
    def __getitem__(self, index):
        if isinstance(index, slice):
            if index.step and index.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self, index.start, index.stop, valueType=self.valueType)
        else:
            return self.itemType(self, index)
    def __len__(self):
        return self.orig.__len__()
    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self._parent, self.brMap, self.itemType)

class CalcVariationsBase:
    """ extra base class for AltCollectionVariations or AltLeafGroupVariations with calculator """
    def __init__(self, withSystName=None):
        self.calc = None
        self.calcProd = None
        self.withSystName = withSystName
    def _initCalc(self, calcProxy, calcHandle=None, args=None):
        self.calc = calcHandle ## handle to the actual module object
        self.calcProd = calcProxy.produce(*args)
    def _initFromCalc(self):
        avail = list(self.calc.available())
        for i,nm in enumerate(avail):
            self.brMapMap[nm] = dict((attN,
                adaptArg(getattr(self.calcProd, attN)(makeConst(i, SizeType))))
                for attN in self.brMapMap[self.withSystName].keys())
        for ky in self.brMapMap[self.withSystName].keys():
            origOp = self.brMapMap[self.withSystName][ky]
            self.brMapMap[self.withSystName][ky] = SystAltOp(
                    self.brMapMap["nominal"][ky], origOp.systName,
                    dict((var, brMap[ky]) for var,brMap in self.brMapMap.items() if var not in (self.withSystName, "nominal", "raw")),
                    valid=[ vr for vr in avail if vr != "nominal" ])

class CalcLeafGroupVariations(AltLeafGroupVariations, CalcVariationsBase):
    """ Set of groups with alternative names for some branches, with calculator """
    def __init__(self, parent, orig, brMapMap, altProxyType, withSystName=None):
        super(CalcLeafGroupVariations, self).__init__(parent, orig, brMapMap, altProxyType)
        CalcVariationsBase.__init__(self, withSystName=withSystName)

class CalcCollectionVariations(AltCollectionVariations, CalcVariationsBase):
    """ Set of collections with alternative names for some branches, with calculator """
    def __init__(self, parent, orig, brMapMap, altItemType=None, withSystName=None):
        super(CalcCollectionVariations, self).__init__(parent, orig, brMapMap, altItemType=altItemType)
        CalcVariationsBase.__init__(self, withSystName=withSystName)

class CombinationProxy(TupleBaseProxy):
    ## NOTE decorated rdfhelpers::Combination
    def __init__(self, cont, parent):
        self.cont = cont
        super(CombinationProxy, self).__init__("struct", parent=parent) ## parent=ObjectProxy for a combination (indices)
    @property
    def _idx(self):
        return self._parent._index
    @property
    def index(self):
        return self._parent.index
    def __getitem__(self, i):
        if isinstance(i, slice):
            if i.step and i.step != 1:
                raise RuntimeError("Slices with non-unit step are not implemented")
            return SliceProxy(self, i.start, i.stop, valueType="struct")
        else:
            idx = makeConst(i, SizeType)
            return self.cont.base(i)[self._parent.result.get(idx)]
    ## TODO add more (maybe simply defer to rdfhelpers::Combination object
    def __repr__(self):
        return "{0}({1!r}, {2!r})".format(self.__class__.__name__, self.cont, self._parent)

class CombinationListProxy(TupleBaseProxy,ListBase):
    ## NOTE list of decorated rdfhelpers::Combination
    ## TODO check if this works with selection (it should...)
    def __init__(self, ranges, parent):
        ListBase.__init__(self) ## FIXME check above how to do this correctly...
        self.ranges = ranges
        super(CombinationListProxy, self).__init__(parent.resultType, parent=parent)
    def _getItem(self, idx):
        return CombinationProxy(self, adaptArg(self._parent.result[idx]))
    def base(self, i):
        return self.ranges[i]._base
    def __len__(self):
        return self._parent.result.size()
    def __repr__(self):
        return "{0}({1!r})".format(self.__class__.__name__, self.ranges, self._parent)
