#pragma once

#include "kinematicvariations.h"

#include "TRandom3.h"
#include <map>
#include <boost/container/flat_map.hpp>
#include "JetResolution.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
class FactorizedJetCorrectorCalculator;

class JetMETVariationsCalculatorBase {
public:
  using p4compv_t = ROOT::VecOps::RVec<float>;

  JetMETVariationsCalculatorBase() = default;
  ~JetMETVariationsCalculatorBase();

  // set up smearing (and JER systematics)
  void setSmearing(const std::string& ptResolution, const std::string& ptResolutionSF, bool doGenMatch, float genMatch_maxDR=-1., float genMatch_maxDPT=-1.)
  {
    m_doSmearing = true;
    m_jetPtRes   = JME::JetResolution(ptResolution);
    m_jetEResSF = JME::JetResolutionScaleFactor(ptResolutionSF);
    m_smearDoGenMatch = doGenMatch;
    m_genMatch_dR2max = genMatch_maxDR*genMatch_maxDR;
    m_genMatch_dPtmax = genMatch_maxDPT;
  }

  void setJEC(const std::vector<JetCorrectorParameters>& jecParams);

  void addJESUncertainty(const std::string& name, const JetCorrectorParameters& params)
  {
    m_jesUncSources.emplace(std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple(params));
  }
protected:
  std::size_t findGenMatch(const float pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const float resolution ) const;

  // config options
  bool m_doSmearing{false}, m_smearDoGenMatch; // default: yes, yes
  float m_genMatch_dR2max, m_genMatch_dPtmax;  // default: R/2 (0.2) and 3
  // parameters and helpers
  JME::JetResolution m_jetPtRes;
  JME::JetResolutionScaleFactor m_jetEResSF;
  mutable TRandom3 m_random; // for resolution
  struct jetcorrdeleter { void operator()(FactorizedJetCorrectorCalculator*) const; };
  // TODO if these would have pure interface functions operator() and produceModifiedCollections could be const (and largely thread-safe)
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrector;
  std::map<std::string,JetCorrectionUncertainty> m_jesUncSources;
  //boost::container::flat_map<std::string,JetCorrectionUncertainty> m_jesUncSources; // problem with
};

class JetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMCollection;
  using LorentzVector = result_t::LorentzVector;

  JetVariationsCalculator() = default;

  std::vector<std::string> availableProducts() const;
  // interface for NanoAOD
  result_t produceModifiedCollections(
      const p4compv_t& jet_rawpt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawfactor, const p4compv_t& jet_area, const float rho,
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
      );
};
