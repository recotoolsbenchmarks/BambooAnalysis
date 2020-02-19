#pragma once

#include "BinnedValuesJSONParser.h"
#include "BinnedValues.h"
#include <stdexcept>

class ILeptonScaleFactor {
public:
  virtual ~ILeptonScaleFactor() {}
  virtual float get(const Parameters& parameters, SystVariation variation) const = 0;
};
class IDiLeptonScaleFactor {
public:
  virtual ~IDiLeptonScaleFactor() {}
  virtual float get(const Parameters& l1Param, const Parameters& l2Param, SystVariation variation) const = 0;
};

class IJetScaleFactor {
public:
  virtual ~IJetScaleFactor() {}

  enum Flavour {
    Light  = 0,
    Charm  = 1,
    Beauty = 2
  };
  static inline Flavour get_flavour(int hadron_flavour) {
    switch(hadron_flavour) {
      case 5:
        return Flavour::Beauty;
      case 4:
        return Flavour::Charm;
      default:
        return Flavour::Light;
    }
  }
  virtual float get(const Parameters& parameters, Flavour flavour, SystVariation variation) const = 0;
};

/**
 * Simple case: one scale factor from file
 */
class ScaleFactor : public ILeptonScaleFactor {
public:
  explicit ScaleFactor(const std::string& file)
  : m_values{BinnedValuesJSONParser::parse_file(file)}
  {}
  virtual ~ScaleFactor() {}

  virtual float get(const Parameters& parameters, SystVariation variation) const override final {
    const auto valForBin = m_values.get(parameters);
    switch (variation) {
      case Nominal:
        return valForBin[Nominal];
      case Up:
        return valForBin[Nominal]+valForBin[Up];
      case Down:
        return valForBin[Nominal]-valForBin[Down];
    }
  }
private:
  BinnedValues m_values;
};

/**
 * Dilepton scalefactor from two lepton scalefactors (takes a parameter set for each)
 */
class DiLeptonFromLegsScaleFactor : public IDiLeptonScaleFactor {
public:
  DiLeptonFromLegsScaleFactor( std::unique_ptr<ILeptonScaleFactor>&& l1_leg1
                             , std::unique_ptr<ILeptonScaleFactor>&& l1_leg2
                             , std::unique_ptr<ILeptonScaleFactor>&& l2_leg1
                             , std::unique_ptr<ILeptonScaleFactor>&& l2_leg2 )
  : m_l1_leg1{std::move(l1_leg1)}
  , m_l1_leg2{std::move(l1_leg2)}
  , m_l2_leg1{std::move(l2_leg1)}
  , m_l2_leg2{std::move(l2_leg2)}
  {}
  virtual ~DiLeptonFromLegsScaleFactor() {}

  virtual float get(const Parameters& l1Param, const Parameters& l2Param, SystVariation variation) const override final {
    const float eff_lep1_leg1 = m_l1_leg1->get(l1Param, variation);
    const float eff_lep1_leg2 = m_l1_leg2->get(l1Param, variation);
    const float eff_lep2_leg1 = m_l2_leg1->get(l2Param, variation);
    const float eff_lep2_leg2 = m_l2_leg2->get(l2Param, variation);

    if ( variation == Nominal ) {
      return -(eff_lep1_leg1 * eff_lep2_leg1) +
          (1 - (1 - eff_lep1_leg2)) * eff_lep2_leg1 +
          eff_lep1_leg1 * (1 - (1 - eff_lep2_leg2)) ;
    } else {
      const float nom_eff_lep1_leg1 = m_l1_leg1->get(l1Param, Nominal);
      const float nom_eff_lep1_leg2 = m_l1_leg2->get(l1Param, Nominal);
      const float nom_eff_lep2_leg1 = m_l2_leg1->get(l2Param, Nominal);
      const float nom_eff_lep2_leg2 = m_l2_leg2->get(l2Param, Nominal);

      const float nominal = -(nom_eff_lep1_leg1 * nom_eff_lep2_leg1) +
          (1 - (1 - nom_eff_lep1_leg2)) * nom_eff_lep2_leg1 +
          nom_eff_lep1_leg1 * (1 - (1 - nom_eff_lep2_leg2)) ;

      const float error_squared =
          std::pow(1 - nom_eff_lep2_leg1 - (1 - nom_eff_lep2_leg2), 2) *
          std::pow(nom_eff_lep1_leg1-eff_lep1_leg1, 2) +
          std::pow(nom_eff_lep2_leg1, 2) *
          std::pow(nom_eff_lep1_leg2-eff_lep1_leg2, 2) +
          std::pow(1 - nom_eff_lep1_leg1 - (1 - nom_eff_lep1_leg2), 2) *
          std::pow(nom_eff_lep2_leg1-eff_lep2_leg1, 2) +
          std::pow(nom_eff_lep1_leg1, 2) *
          std::pow(nom_eff_lep2_leg2-eff_lep2_leg2, 2);

      if ( variation == Up ) {
        return std::min(nominal+std::sqrt(error_squared), 1.f);
      } else if ( variation == Down ) {
        return std::max(0.f, nominal-std::sqrt(error_squared));
      } else {
        throw std::invalid_argument("Unsupported variation: "+std::to_string(variation));
      }
    }
  }
private:
  std::unique_ptr<ILeptonScaleFactor> m_l1_leg1;
  std::unique_ptr<ILeptonScaleFactor> m_l1_leg2;
  std::unique_ptr<ILeptonScaleFactor> m_l2_leg1;
  std::unique_ptr<ILeptonScaleFactor> m_l2_leg2;
};

/**
 * B-tagging scale factor (depends on flavour)
 */
class BTaggingScaleFactor : public IJetScaleFactor {
public:
  BTaggingScaleFactor(const std::string& lightFile, const std::string& charmFile, const std::string& beautyFile)
  : m_lightValues {BinnedValuesJSONParser::parse_file(lightFile )}
  , m_charmValues {BinnedValuesJSONParser::parse_file(charmFile )}
  , m_beautyValues{BinnedValuesJSONParser::parse_file(beautyFile)}
  {}
  virtual ~BTaggingScaleFactor() {}

  virtual float get(const Parameters& parameters, Flavour flavour, SystVariation variation) const override final {

    std::vector<float> values;
    switch(flavour) {
      case Flavour::Beauty:
        values = m_beautyValues.get(parameters);
        break;
      case Flavour::Charm:
        values = m_charmValues.get(parameters);
        break;
      case Flavour::Light:
        values = m_lightValues.get(parameters);
        break;
      default:
        throw std::invalid_argument("Unsupported flavour: "+std::to_string(flavour));
    }
    switch(variation) {
      case Nominal:
        return values[Nominal];
      case Up:
        return values[Nominal]+values[Up];
      case Down:
        return values[Nominal]-values[Down];
    }
  }
private:
  BinnedValues m_lightValues;
  BinnedValues m_charmValues;
  BinnedValues m_beautyValues;
};

#include <vector>
#include <algorithm>
#include "range.h"
/**
 * Weight between different scale factors
 */
template<class ISingle,typename... Args>
class WeightedScaleFactor : public ISingle {
public:
  WeightedScaleFactor( const std::vector<float>& probs, std::vector<std::unique_ptr<ISingle>>&& sfs )
  : m_terms{std::move(sfs)}
  {
    const double norm{1./std::accumulate(std::begin(probs), std::end(probs), 0.)};
    std::transform(std::begin(probs), std::end(probs), std::back_inserter(m_probs), [norm] ( float p ) -> float { return norm*p; } );
  }
  virtual ~WeightedScaleFactor() {}

  virtual float get(Args&&... args, SystVariation variation) const override final {
    const auto myRange = rdfhelpers::IndexRange<std::size_t>(m_terms.size());
    return std::accumulate(myRange.begin(), myRange.end(), 0.,
        [this,&args...,variation] ( float wsum, std::size_t i ) -> float {
          return wsum + m_probs[i]*m_terms[i]->get(std::forward<Args>(args)..., variation);
        }
      );
  }
private:
  std::vector<float> m_probs;
  std::vector<std::unique_ptr<ISingle>> m_terms;
};
using WScaleFactor         = WeightedScaleFactor<ILeptonScaleFactor,const Parameters&>;
using WLLScaleFactor       = WeightedScaleFactor<IDiLeptonScaleFactor,const Parameters&,const Parameters&>;
using WBTaggingScaleFactor = WeightedScaleFactor<IJetScaleFactor,const Parameters&,IJetScaleFactor::Flavour>;

#include "bamboorandom.h"
/**
 * Sample according to fractions from different scale factors
 */
template<class ISingle,typename... Args>
class SampledScaleFactor : public ISingle
{
public:
  SampledScaleFactor( const std::vector<float> probs, std::vector<std::unique_ptr<ISingle>>&& sfs )
  : m_probs{probs}
  , m_terms{std::move(sfs)}
  {}
  virtual ~SampledScaleFactor() {}

  virtual float get(Args... args, std::uint32_t seed, SystVariation variation) const override final {
    auto& rg = rdfhelpers::getStdMT19937(seed);
    return m_terms[std::discrete_distribution<std::size_t>(m_probs)(rg)]->get(std::forward<Args>(args)..., variation);
  }
private:
  std::vector<float> m_probs;
  std::vector<std::unique_ptr<ISingle>> m_terms;
};
using SmpScaleFactor         = SampledScaleFactor<ILeptonScaleFactor,const Parameters&>;
using SmpLLScaleFactor       = SampledScaleFactor<IDiLeptonScaleFactor,const Parameters&,const Parameters&>;
using SmpBTaggingScaleFactor = SampledScaleFactor<IJetScaleFactor,const Parameters&,IJetScaleFactor::Flavour>;
