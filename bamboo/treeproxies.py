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

boolType = "bool"
NumberTypes = set(("Float_t", "Double_t", "Int_t", "UInt_t", "Bool_t", "Char_t", "UChar_t", "ULong64_t", "int", "unsigned", "unsigned short", "char", "signed char", "unsigned char", "float", "double", "Short_t", "size_t", "std::size_t", "unsigned short", "bool", "unsigned long")) ## there are many more (at least unsigned)
import re
vecPat = re.compile("(?:vector|ROOT\:\:VecOps\:\:RVec)\<(?P<item>[a-zA-Z_0-9\<\>,\: ]+)\>")

def makeProxy(typeName, parent, length=None):
    if length is not None:
        return ArrayProxy(parent, typeName, length)
    if typeName in NumberTypes:
        return NumberProxy(parent, typeName)
    else:
        m = vecPat.match(typeName)
        if m:
            return VectorProxy(parent, typeName)
        else:
            return ObjectProxy(parent, typeName)
def makeConst(value, typeHint):
    return makeProxy(typeHint, adaptArg(value, typeHint))

class NumberProxy(TupleBaseProxy):
    """ Proxy for a number (integer or floating-point) """
    def __init__(self, parent, typeName):
        super(NumberProxy, self).__init__(typeName, parent=parent)
    def __repr__(self):
        return "NumberProxy({0!r}, {1!r})".format(self._parent, self._typeName)
    def _binaryOp(self, opName, other, outType="Double_t"):
        return MathOp(opName, self, other, outType=outType).result
## operator overloads
for nm,opNm in {
          "__add__" : "add"
        , "__sub__" : "subtract"
        , "__mul__" : "multiply"
        , "__truediv__" : " divide"
        , "__div__" : "divide"
        , "__nonzero__" : "nonzero"
        }.items():
    setattr(NumberProxy, nm, (lambda oN : (lambda self, other : self._binaryOp(oN, other)))(opNm))
for nm in ("__lt__", "__le__", "__eq__", "__ne__", "__gt__", "__ge__"):
    setattr(NumberProxy, nm, (lambda oN, oT : (lambda self, other : self._binaryOp(oN, other, outType=oT)))(nm.strip("_"), boolType))

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
    def deps(self, defCache=None, select=(lambda x : True)):
        yield from []

class ListBase(object):
    """ Interface definition for range proxies (Array/Vector, split object vector, selection/reordering) """
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        yield from self._size.deps(defCache=defCache, select=select)
    def __eq__(self, other):
        return isinstance(other, ContainerGroupProxy) and ( self._size == other._size ) and ( self.valuetype == other.valuetype ) and ( self._parent == other._parent ) and ( self._prefix == other._prefix )

class ObjectProxy(NumberProxy):
    """
    Imitate an object
    """
    __slots__ = ("_typ",)
    def __init__(self, parent, typeName):
        import ROOT
        self._typ = getattr(ROOT, typeName) ## NOTE could also use TClass machinery
        super(ObjectProxy, self).__init__(parent, typeName)
    def __getattr__(self, name):
        if name not in dir(self._typ):
            raise AttributeError("Type {0} has no member {1}".format(self._typeName, name))
        import ROOT
        if hasattr(self._typ, name) and isinstance(getattr(self._typ, name), ROOT.MethodProxy):
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
        ## TODO maybe tihs is a good place to resolve the right overload? or do some arguments checking
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
        ## TODO maybe tihs is a good place to resolve the right overload? or do some arguments checking
        return CallMethod(self._name, tuple(args)).result
    def __repr__(self):
        return "MethodProxy({0!r})".format(self._name)

class VectorProxy(ObjectProxy,ListBase):
    """ Vector-as-array (to be eliminated with var-array branches / generalised into object) """
    def __init__(self, parent, typeName):
        ListBase.__init__(self)
        import ROOT
        vecClass = getattr(ROOT, typeName)
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
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        yield from self.parent.deps(defCache=defCache, select=select)
    def __eq__(self, other):
        return isinstance(other, VectorProxy) and ( self._parent == other._parent ) and ( self.valuetype == other.valuetype )

class SelectionProxy(TupleBaseProxy,ListBase):
    """ Proxy for a selection from an iterable (ContainerGroup/ other selection etc.) """
    def __init__(self, parent, base):
        ListBase.__init__(self)
        self.valueType = base.valueType
        self._base = base
        ## the list of indices is stored as the parent
        super(SelectionProxy, self).__init__(base.valueType, parent=parent)
    @property
    def _idxs(self):
        return self._parent.result
    def __getitem__(self, index):
        return self._base[GetItem(self._idxs, SizeType, index)]
    def __len__(self):
        return self._idxs.__len__()
    def __repr__(self):
        return "SelectionProxy({0!r}, {1!r})".format(self._parent, self._base)
    #
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (self._parent, self._base):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
    def __eq__(self, other):
        return isinstance(other, SelectionProxy) and ( self._parent == other._parent ) and ( self.valuetype == other.valuetype ) and ( self._base == other._base )

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
        return self.op[inst._parent._parent.indices()[inst._idx]]
class varItemRefProxy(object):
    def __init__(self, op, getTarget):
        self.op = op
        self.getTarget = getTarget
    def __get__(self, inst, cls):
        return self.getTarget(inst)[self.op[inst._parent._parent.indices()[inst._idx]]]

class Variations(TupleBaseProxy):
    """ Proxy with variations (that can be dynamically added) """
    def __init__(self, parent, nominal, varItemType=None):
        self._vars = dict()
        self._vars["nominal"] = nominal
        self._varItemType = varItemType if varItemType else nominal.valuetype
        super(Variations, self).__init__("dict", parent=parent)
    @property
    def nominal(self):
        return self._vars["nominal"]
    def __getitem__(self, key):
        return self._vars[key]
    def register(self, key, modified):
        if key in self._vars:
            raise ValueError("Variation '{0}' is already registered".format(key))
        self._vars[key] = modified

class ModifiedCollectionProxy(TupleBaseProxy,ListBase):
    """ Collection with a selection and modified branches """
    def __init__(self, parent, base, itemType=None):
        ## parent has a member indices() and one for each of modBranches
        ListBase.__init__(self)
        self.valueType = base.valueType ## for ListBase
        self._base = base
        self.itemType = itemType ## for actual items
        super(ModifiedCollectionProxy, self).__init__(self.valueType, parent=parent)
    @property
    def _idxs(self):
        return self._parent.indices()
    def __len__(self):
        return self._idxs.__len__()
    def __getitem__(self, index):
        return self.itemType(self, index)
    def __repr__(self):
        return "{0}({1!r}, {2!r}, {3!r})".format(self.__class__.__name__, self._parent, self._base, self.itemType)
    ##
    def deps(self, defCache=cppNoRedir, select=(lambda x : True)):
        for arg in (self._parent, self._base):
            if select(arg):
                yield arg
            yield from arg.deps(defCache=defCache, select=select)
    def __eq__(self, other):
        return isinstance(other, self.__class__) and ( self._parent == other._parent ) and ( self._base == other._base ) and ( self.itemType == other.itemType )
