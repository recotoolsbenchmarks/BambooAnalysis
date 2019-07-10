"""
The :py:mod:`bamboo.scalefactors` module contains helper methods
for configuring scale factors, fake rates etc.

The basic configuration parameter is the JSON file path for a set of scalefactors.
There two basic types are

- lepton scale factors (dependent on a number of object variables, e.g. pt and eta),
- jet (b-tagging) scale factors (grouped set for different flavours, for convenience)

Different values (depending on the data-taking period)
can be taken into account by weighting or by randomly sampling.
"""
__all__ = ("get_scalefactor", "lumiPerPeriod")

from itertools import chain

from . import treefunctions as op

#: Integrated luminosity (in 1/pb) per data taking period
lumiPerPeriod_default = {
      "Run2016B" : 5785.152 ## averaged 5783.740 (DoubleMuon), 5787.976 (DoubleEG) and 5783.740 (MuonEG) - max dev. from average is 5.e-4 << lumi syst
    , "Run2016C" : 2573.399
    , "Run2016D" : 4248.384
    , "Run2016E" : 4009.132
    , "Run2016F" : 3101.618
    , "Run2016G" : 7540.488
    , "Run2016H" : 8605.689
    ##
    , "Run271036to275783" : 6274.191
    , "Run275784to276500" : 3426.131
    , "Run276501to276811" : 3191.207
    }


# TODO maybe move this elsewhere
binningVariables_nano = {
      "Eta" : lambda obj : obj.p4.Eta()
    , "AbsEta" : lambda obj : op.abs(obj.p4.Eta())
    , "ClusEta" : lambda el : el.eta+el.deltaEtaSC
    , "AbsClusEta" : lambda el : op.abs(el.eta+el.deltaEtaSC)
    , "Pt" : lambda obj : obj.p4.Pt()
    }

def getBinningVarNames(jsonpath):
    import json
    with open(jsonpath, "r") as jsf:
        cont = json.load(jsf)
    return tuple(cont["variables"])

class BinningParameters(object):
    def __init__(self, binningVars):
        self.binningVars = binningVars
    def __call__(self, obj):
        return op.construct("Parameters",
                       (op.initList("std::initializer_list<Parameters::value_type::value_type>", "Parameters::value_type::value_type", (
                           op.initList("Parameters::value_type::value_type", "float", (op.extVar("int", "BinningVariable::{0}".format(bvNm.replace("ClusEta", "Eta"))), bv(obj)))
                           for bvNm,bv in self.binningVars.items())),)
                   )

def getBinningParameters(bVarNames, isElectron=False, moreVars=dict(), paramDefs=dict()):
    if isElectron:
        bVarNames = [ k.replace("Eta", "ClusEta") for k in bVarNames ]
    theDict = dict(paramDefs)
    theDict.update(moreVars)
    return BinningParameters(dict((k,theDict[k]) for k in bVarNames))

class ScaleFactor(object):
    def __init__(self, cppDef=None, args=None, iface="ILeptonScaleFactor"):
        self._cppDef = cppDef
        self._args = args
        self.sfOp = op.define(iface, cppDef)
    def __call__(self, obj, variation="Nominal"):
        from .treedecorators import makeConst, boolType
        expr = self.sfOp.get(*tuple(chain(
                   list(a(obj) for a in self._args)
                 , (op.extVar("int", variation),)
               )))
        return expr

def get_scalefactor(objType, key, periods=None, combine=None, additionalVariables=dict(), sfLib=dict(), paramDefs=dict(), lumiPerPeriod=lumiPerPeriod_default):
    """ Construct a scalefactor callable

    :param objType: object type: ``"lepton"``, ``"dilepton"``, or ``"jet"``
    :param key: key in ``sfLib`` (or tuple of keys, in case of a nested dictionary)
    :param periods: data taking periods to consider when combining different scalefactors
    :param combine: combination strategy (``"weight"`` or ``"sample"``)
    :param paramDefs: dictionary of binning variable definitions (name to callable)
    :param additionalVariables: additional binning variable definitions (TODO: remove)
    :param lumiPerPeriod: alternative definitions and relative weights of run periods

    :returns: a callable that takes ``(object, variation="Nominal")`` and returns a floating-point number proxy
    """
    ##
    ## Interpret args, get defaults etc.
    ## 
    if isinstance(key, tuple):
        # interpret key=("a", "b") as ...["a"]["b"]
        mainKey = key[0]
        config = sfLib[key[0]]
        for idx in range(1,len(key)):
            config = config[key[idx]]
    else:
        mainKey = key
        config = sfLib[key]

    if periods is None:
        if "2016" in mainKey:
            periods = [ "Run2016{0}".format(ltr) for ltr in "BCDEFGH" ]
        else: ## TODO similar for 2017, 2018
            periods = []
    periods = set(periods)
    
    if combine is not None:
        combPrefix = { "weight" : "W"
                     , "sample" : "Smp" }.get(combine, "W")

    ##
    ## Construct scalefactors
    ##
    if objType == "lepton":
        iface = "ILeptonScaleFactor"
        isElectron = (key[0].split("_")[0] == "electron")
        if isinstance(config, str):
            return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(config),
                    args=(getBinningParameters(getBinningVarNames(config), isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                    iface=iface)
        else:
            if combPrefix == "":
                raise ValueError("A combination mode needs to be specified for this scale factor")
            selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),path)
                    for eras,path in config if any(ier in periods for ier in eras))))
            if len(selConfigs) < 1:
                raise RuntimeError("Zero period configs selected for config {0} with periods {1}".format(", ".join("({0} : {1})".format(list(eras), path) for eras, path in config), list(periods)))
            elif len(selConfigs) == 1:
                return ScaleFactor(cppDef='const ScaleFactor <<name>>{{"{0}"}};'.format(selConfigs[0][1]),
                        args=(getBinningParameters(getBinningVarNames(selConfigs[0][1]), isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                        iface=iface)
            else:
                bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,iPth in selConfigs))
                return ScaleFactor(cppDef=(
                            'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                'std::make_unique<ScaleFactor>("{0}")'.format(path) for wgt, path in selConfigs), iface=iface)+
                            'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,path in selConfigs), cmb=combPrefix)+
                              'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                            ),
                        args=(getBinningParameters(bVarNames, isElectron=isElectron, moreVars=additionalVariables, paramDefs=paramDefs),),
                        iface=iface)
    elif objType == "dilepton":
        iface = "IDiLeptonScaleFactor"
        if isinstance(config, tuple) and len(config) == 4:
            if not all(isinstance(iCfg, str) for iCfg in config):
                raise TypeError("Config for dilepton scale factor should be quadruplet of paths or list f weights and triplets, found {0}".format(config))

            return ScaleFactor(cppDef="const DiLeptonFromLegsScaleFactor <<name>>{{{0}}};".format(", ".join(
                        'std::make_unique<ScaleFactor>("{0}")'.format(leplepCfg) for leplepCfg in config)),
                    args=[ (lambda bp : (lambda ll : bp(ll[0])))(getBinningParameters(set(chain(getBinningVarNames(config[0]), getBinningVarNames(config[1]))), moreVars=additionalVariables, paramDefs=paramDefs))
                         , (lambda bp : (lambda ll : bp(ll[1])))(getBinningParameters(set(chain(getBinningVarNames(config[2]), getBinningVarNames(config[3]))), moreVars=additionalVariables, paramDefs=paramDefs)) ],
                    iface=iface)
        else:
            raise NotImplementedError("Still to do this part")
    elif objType == "jet":
        iface = "IJetScaleFactor"
        if isinstance(config, tuple) and len(config) == 3:
            if not all(isinstance(iCfg, str) for iCfg in config):
                raise TypeError("Config for b-tagging should be triplet of paths or list of weights and triplets, found {0}".format(config))
            else:
                bVarNames = set(chain.from_iterable(getBinningVarNames(iCfg) for iCfg in config))
                return ScaleFactor(cppDef='const BTaggingScaleFactor <<name>>{{{0}}};'.format(", ".join('"{0}"'.format(iCfg) for iCfg in config)),
                        args=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                        iface=iface)
        else:
            if not ( all((isinstance(iCfg, tuple) and len(iCfg) == 3 and all(isinstance(iPth, str) for iPth in iCfg) ) for iCfg in config) ):
                raise TypeError("Config for b-tagging should be triplet of paths or list of weights and triplets, found {0}".format(config))
            else:
                if combPrefix == "":
                    raise ValueError("A combination mode needs to be specified for this scale factor")
                selConfigs = list(filter((lambda elm : elm[0] != 0.), # only keep those with nonzero lumi
                    ((sum(lumiPerPeriod[ier] for ier in eras if ier in periods),paths)
                        for eras,paths in config if any(ier in periods for ier in eras))))
                if len(selConfigs) < 1:
                    raise RuntimeError("Zero configs")
                elif len(selConfigs) == 1:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iCfg) for iCfg in selConfigs[0]))
                    return ScaleFactor(cppDef='const BTaggingScaleFactor <<name>>{{{0}}};'.format(", ".join('"{0}"'.format(iCfg) for iCfg in selConfigs[0])),
                            args=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                            iface=iface)
                else:
                    bVarNames = set(chain.from_iterable(getBinningVarNames(iPth) for iWgt,paths in selConfigs for iPth in paths))
                    return ScaleFactor(cppDef=(
                                'std::unique_ptr<{iface}> tmpSFs_<<name>>[] = {{ {0} }};\n'.format(", ".join(
                                    'std::make_unique<BTaggingScaleFactor>({0})'.format(", ".join('"{0}"'.format(iPth) for iPth in paths)) for wgt, paths in selConfigs), iface=iface)+
                                'const {cmb}ScaleFactor <<name>>{{ {{ {0} }}, '.format(", ".join("{0:e}".format(wgt) for wgt,paths in selConfigs), cmb=combPrefix)+
                                  'std::vector<std::unique_ptr<{iface}>>{{std::make_move_iterator(std::begin(tmpSFs_<<name>>)), std::make_move_iterator(std::end(tmpSFs_<<name>>))}} }};'.format(iface=iface)
                                ),
                            arg=(getBinningParameters(bVarNames, moreVars=additionalVariables, paramDefs=paramDefs), (lambda j : op.extMethod("IJetScaleFactor::get_flavour")(j.hadronFlavor))),
                            iface=iface)
    else:
        raise ValueError("Unknown object type: {0}".format(objType))
