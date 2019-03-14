#pragma once

#include <ROOT/RVec.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>

namespace rdfhelpers {

class ModifiedKinCollection {
public:
  using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
  using Indices = ROOT::VecOps::RVec<std::size_t>;
  using Momenta  = ROOT::VecOps::RVec<LorentzVector>;

  ModifiedKinCollection() = default;
  ModifiedKinCollection( const Indices& indices, const Momenta& momenta )
    : m_indices(indices), m_momenta(momenta) {}
  ModifiedKinCollection( Indices&& indices, Momenta&& momenta )
    : m_indices(indices), m_momenta(momenta) {}
  const Indices& indices() const { return m_indices; }
  const Momenta& momenta() const { return m_momenta; }
private:
  Indices m_indices;
  Momenta m_momenta;
}; // TODO make it jet-specific and add MET

};

#include <random>
#include <map>
#include <boost/container/flat_map.hpp>
#include "JetResolution.h"
#include "JetCorrectorParameters.h"
#include "JetCorrectionUncertainty.h"
class FactorizedJetCorrectorCalculator;

class JMESystematicsCalculator {
public:
  using result_entry_t = rdfhelpers::ModifiedKinCollection;
  using result_t = boost::container::flat_map<std::string,result_entry_t>;
  using p4compv_t = ROOT::VecOps::RVec<float>;
  using LorentzVector = result_entry_t::LorentzVector;

  JMESystematicsCalculator() {}

  ~JMESystematicsCalculator();

  // set up smearing (and JER systematics)
  void setSmearing(const std::string& ptResolution, const std::string& ptResolutionSF, bool doGenMatch, float genMatch_maxDR=-1., float genMatch_maxDPT=-1.)
  {
    m_doSmearing = true;
    m_jetPtRes   = JME::JetResolution(ptResolution);
    m_jetEResSF = JME::JetResolutionScaleFactor(ptResolutionSF);
    m_smearDoGenMatch = doGenMatch;
    m_genMatch_dRmax = genMatch_maxDR;
    m_genMatch_dPtmax = genMatch_maxDPT;
  }

  void setJEC(const std::vector<JetCorrectorParameters>& jecParams);

  void addJESUncertainty(const std::string& name, const JetCorrectorParameters& params)
  {
    m_jesUncSources.emplace(std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple(params));
  }
  // TODO retire in favour of the one above?
  void addJESUncertainty(const std::string& name, const std::string& fileName, const std::string& section="__default__")
  {
    m_jesUncSources.emplace(std::piecewise_construct,
        std::forward_as_tuple(name),
        std::forward_as_tuple(JetCorrectorParameters(fileName, ( section != "__default__" ? section : name ))));
  }

  bool hasProduct(const std::string& name) const;

  // interface for NanoAOD
  result_t produceModifiedCollections(
      const p4compv_t& jet_rawpt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
      const p4compv_t& jet_rawfactor, const p4compv_t& jet_area, const float rho,
      const float met_phi, const float met_pt, const float sumEt,
      // TODO unclustered?
      const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
      );
  // TODO interface for Framework+TTWAnalysis
  // result_t operator() ( const std::vector<JMESystematicsCalculator::LorentzVector>& jet_p4, const std::vector<float>& jet_rawfactor, const std::vector<float>& jet_area, const float rho, const std::vector<JMESystematicsCalculator::LorentzVector>& genjet_p4 ) const;

private:
  struct Jet {
    std::size_t i;
    LorentzVector p4;
    Jet(std::size_t i_, const LorentzVector& p4_) : i(i_), p4(p4_) {}
    Jet(std::size_t i_, float pt, float eta, float phi, float m) : i(i_), p4(pt, eta, phi, m) {}
  };
  void sort(std::vector<Jet>&) const;
  result_entry_t convertToModifKin( const std::vector<Jet>& jets ) const;
  std::size_t findGenMatch( const LorentzVector& p4, const ROOT::VecOps::RVec<LorentzVector>& genp4, const float resolution ) const;

  result_t produceModifiedCollections( const ROOT::VecOps::RVec<LorentzVector>& jetp4, const ROOT::VecOps::RVec<double>& jet_rawfactor, const p4compv_t& jet_area, const float rho, const ROOT::VecOps::RVec<LorentzVector>& genp4 );

  // config options
  bool m_doSmearing{false}, m_smearDoGenMatch; // TODO default: yes, yes
  float m_genMatch_dRmax, m_genMatch_dPtmax;   // TODO default: R/2 (0.2) and 3
  // parameters and helpers
  JME::JetResolution m_jetPtRes;
  JME::JetResolutionScaleFactor m_jetEResSF;
  mutable std::mt19937 m_random; // for resolution
  struct jetcorrdeleter { void operator()(FactorizedJetCorrectorCalculator*) const; };
  // TODO if these would have pure interface functions operator() and produceModifiedCollections could be const (and largely thread-safe)
  std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter> m_jetCorrector;
  std::map<std::string,JetCorrectionUncertainty> m_jesUncSources;
  //boost::container::flat_map<std::string,JetCorrectionUncertainty> m_jesUncSources; // problem with
};
