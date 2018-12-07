#pragma once

#include <ROOT/RVec.hxx>
#include <Math/LorentzVector.h>

namespace rdfhelpers {

class ModifiedKinCollection {
public:
  using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
  using Indices = ROOT::VecOps::RVec<std::size_t>;
  using Momenta  = ROOT::VecOps::RVec<LorentzVector>;

  ModifiedKinCollection( const Indices& indices, const Momenta& momenta )
    : m_indices(indices), m_momenta(momenta) {}
  ModifiedKinCollection( Indices&& indices, Momenta&& momenta )
    : m_indices(indices), m_momenta(momenta) {}
  const Indices& indices() const { return m_indices; }
  const Momenta& momenta() const { return m_momenta; }
private:
  Indices m_indices;
  Momenta m_momenta ;
};

template<typename RANGE, typename FUN, typename PRED>
ModifiedKinCollection modifyKinCollection(const RANGE& range, FUN&& fun, PRED&& pred)
{
  ModifiedKinCollection::Indices indices;
  ModifiedKinCollection::Momenta momenta;
  for ( const auto obj : range ) {
    ModifiedKinCollection::LorentzVector mod = fun(obj);
    if ( pred(mod, obj) ) {
      indices.push_back(typename RANGE::value_type(obj));
      momenta.push_back(ModifiedKinCollection::LorentzVector(mod));
    }
  }
  return ModifiedKinCollection(std::move(indices), std::move(momenta));
}

};
