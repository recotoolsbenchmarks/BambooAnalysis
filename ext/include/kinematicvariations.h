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
}; // TODO make a jet-specific version with MET

};
