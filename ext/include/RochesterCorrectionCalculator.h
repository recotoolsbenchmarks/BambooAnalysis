#pragma once

#include "kinematicvariations.h"
#include <boost/container/flat_map.hpp>
#include <memory>
#include <random>
#include <Rtypes.h>
class RoccoR;

class RochesterCorrectionCalculator {
public:
  using result_entry_t = rdfhelpers::ModifiedKinCollection;
  using result_t = boost::container::flat_map<std::string,result_entry_t>;
  using p4compv_t = ROOT::VecOps::RVec<float>;
  using LorentzVector = result_entry_t::LorentzVector;

  RochesterCorrectionCalculator() {}

  ~RochesterCorrectionCalculator();

  void setRochesterCorrection(const std::string& params);

  bool hasProduct(const std::string& name) const;

  // interface for NanoAOD
  result_t produceModifiedCollections(
      const p4compv_t& muon_pt, const p4compv_t& muon_eta, const p4compv_t& muon_phi, const p4compv_t& muon_mass,
      const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const ROOT::VecOps::RVec<Int_t>& muon_genIdx, const p4compv_t& gen_pt) const;

  struct Muon {
    std::size_t i;
    LorentzVector p4;
    Muon(std::size_t i_, const LorentzVector& p4_) : i(i_), p4(p4_) {}
    Muon(std::size_t i_, float pt, float eta, float phi, float m) : i(i_), p4(pt, eta, phi, m) {}
  };

private:
  mutable std::mt19937 m_random; // for resolution
  struct roccordeleter { void operator()(RoccoR*) const; };
  std::unique_ptr<RoccoR,roccordeleter> m_roccor;

  result_t produceModifiedCollections(std::vector<Muon>&& muons, const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const p4compv_t& muon_genpt) const;
};
