#pragma once

#include "Math/LorentzVector.h"
#include "Math/PtEtaPhiE4D.h"

namespace rdfhelpers {

template<typename N1,typename N2,typename N3>
bool in_range( N1 lower, N2 value, N3 upper )
{
  return ( lower < value ) && ( value < upper );
}

template<typename LorentzVector>
ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > withMass( const LorentzVector& p4, float mass )
{
  return ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> >{p4.Pt(), p4.Eta(), p4.Phi(), mass};
}

};
