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
  // TODO if these would have pure interface functions operator() and produce could be const (and largely thread-safe)
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrector;
  std::map<std::string,JetCorrectionUncertainty> m_jesUncSources;
  //boost::container::flat_map<std::string,JetCorrectionUncertainty> m_jesUncSources; // problem with
};

class JetVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedPtMCollection;

  JetVariationsCalculator() = default;

  std::vector<std::string> available() const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
      );
};

class Type1METVariationsCalculator : public JetMETVariationsCalculatorBase {
public:
  using result_t = rdfhelpers::ModifiedMET;

  Type1METVariationsCalculator() = default;

  // additional settings: L1-only JEC
  void setL1JEC(const std::vector<JetCorrectorParameters>& jecParams);
  void setUnclusteredEnergyTreshold(float threshold) { m_unclEnThreshold = threshold; }

  std::vector<std::string> available() const;
  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt,
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
      );
protected:
  float m_unclEnThreshold = 15.;
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1;
  void addVariations(Type1METVariationsCalculator::result_t& out,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<bool>& jet_mask,
      const float rho,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass);
};

class FixEE2017Type1METVariationsCalculator : public Type1METVariationsCalculator {
public:
  FixEE2017Type1METVariationsCalculator() = default;

  // additional settings: full and L1-only JEC used in production
  void setJECProd(const std::vector<JetCorrectorParameters>& jecParams);
  void setL1JECProd(const std::vector<JetCorrectorParameters>& jecParams);

  // interface for NanoAOD
  result_t produce(
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
      const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
      const float rho,
      // MC-only
      const std::uint32_t seed,
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
      // MET-specific
      const float rawmet_phi, const float rawmet_pt, // "RawMET"
      const float met_unclustenupdx, const float met_unclustenupdy,
      const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
      const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF,
      const float defmet_phi, const float defmet_pt, // "MET"
      const float t1met_phi, const float t1met_pt    // "METFixEE2017"
      );
protected:
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorProd;
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrectorL1Prod;
  std::array<float,4> calculateFixEE2017Offset(ROOT::VecOps::RVec<bool>& jet_mask,
      const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
      const float rho);
};
