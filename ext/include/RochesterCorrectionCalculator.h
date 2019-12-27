#pragma once

#include "kinematicvariations.h"
#include <memory>
#include <random>
#include <Rtypes.h>
class RoccoR;

class RochesterCorrectionCalculator {
public:
  using result_t = rdfhelpers::ModifiedPtCollection;
  using p4compv_t = result_t::compv_t;

  RochesterCorrectionCalculator() {}

  ~RochesterCorrectionCalculator();

  void setRochesterCorrection(const std::string& params);

  std::vector<std::string> available() const;

  // interface for NanoAOD
  result_t produce(
      const p4compv_t& muon_pt, const p4compv_t& muon_eta, const p4compv_t& muon_phi, const p4compv_t& muon_mass,
      const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const ROOT::VecOps::RVec<Int_t>& muon_genIdx, const p4compv_t& gen_pt) const;
private:
  mutable std::mt19937 m_random; // for resolution
  struct roccordeleter { void operator()(RoccoR*) const; };
  std::unique_ptr<RoccoR,roccordeleter> m_roccor;
};
