#pragma once

#include "kinematicvariations.h"
#include <boost/container/flat_map.hpp>
#include <memory>
#include <random>
#include <Rtypes.h>
class RoccoR;

class RochesterCorrectionCalculator {
public:
  using result_t = rdfhelpers::ModifiedPtCollection;
  using LorentzVector = result_t::LorentzVector;
  using p4compv_t = result_t::p4compv_t;

  RochesterCorrectionCalculator() {}

  ~RochesterCorrectionCalculator();

  void setRochesterCorrection(const std::string& params);

  std::vector<std::string> availableProducts() const;

  // interface for NanoAOD
  result_t produceModifiedCollections(
      const p4compv_t& muon_pt, const p4compv_t& muon_eta, const p4compv_t& muon_phi, const p4compv_t& muon_mass,
      const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const ROOT::VecOps::RVec<Int_t>& muon_genIdx, const p4compv_t& gen_pt) const;
private:
  mutable std::mt19937 m_random; // for resolution
  struct roccordeleter { void operator()(RoccoR*) const; };
  std::unique_ptr<RoccoR,roccordeleter> m_roccor;
};
