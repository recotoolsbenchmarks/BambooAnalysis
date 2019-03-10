#include "JMESystematicsCalculator.h"

#include "FactorizedJetCorrector.h"

JMESystematicsCalculator::setJEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrector = std::make_unique<FactorizedJetCorrector>(jecParams);
  }
}

JMESystematicsCalculator::jetcorrdeleter::operator() (FactorizedJetCorrector* ptr) const { delete ptr; }

JMESystematicsCalculator::~JMESystematicsCalculator() {}

JMESystematicsCalculator::result_t JMESystematicsCalculator::operator() const (
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const rho_t rho,
    const p4comp_t met_phi, const p4comp_t met_pt, const p4comp_t sumEt,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass )
{
  ROOT::RVec<LorentzVector> jetP4;
  jets.reserve(jet_pt.size());
  for ( std::size_t i{0}; i < jet_pt.size(); ++i ) {
    jets.emplace_back(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
  }
  ROOT::RVec<LorentzVector> genP4;
  genP4.reserve(genjet_pt.size());
  for ( std::size_t i{0}; i < genjet_pt.size(); ++i ) {
    genP4.emplace_back(genjet_pt[i], genjet_eta[i], genjet_phi[i], genjet_mass[i]);
  }
  return produceModifiedCollections(jets, rho, genP4);
}

namespace {
  Float jetSmearFactor( Float_t pt, Float eOrig, Float_t genPt, Float_t ptRes, Float_t sfUncert, std::mt19937& randgen )
  {
    Float_t smear = 1.
    if ( genPt > 0. ) {
      smear = 1. + (sfUncert-1.)*(pt - genPt)/pt;
    } else if ( sfUncert > 1. ) {
      std::normal_distribution<> d(0, ptRes*std::sqrt(sfUncert*sfUncert-1.));
      smear = 1.+d(randgen);
    }
    if ( smear*eOrig < 1.e-2 ) {
      smear = 1.e-2/eOrig;
    }
    return smear;
  }

  void sort( std::vector<Jet>& jets )
  {
    std::sort(std::begin(jets), std::end(jets), [] ( const Jet& j1, const Jet& j2 ) { return j1.Pt() > j2.Pt(); });
  }

  // TODO with orig MET and jets (sumpx,sumpy): calc modif MET(sig), produce bigger results type
  // TODO non-const ref and include sort then
  rdfhelpers::ModifiedKinCollection convertToModifKin( const std::vector<Jet>& jets )
  {
    ModifiedKinCollection::Indices idx;
    ModifiedKinCollection::Momenta mom;
    idx.reserve(jets.size());
    mom.reserve(jets.size());
    for ( const auto& j : jets ) {
      idx.push_back(j.i);
      mom.push_back(j.p4);
    }
    return ModifiedKinCollection(std::move(idx), std::move(mom));
  }
}

std::size_t JMESystematicsCalculator::findGenMatch( const JMESystematicsCalculator::LorentzVector& p4, const ROOT::RVec<JMESystematicsCalculator::LorentzVector>& gen_p4, Float_t resolution )
{
  auto drMin = std::numeric_limits<p4comp_t>::max();
  std::size_t igBest{gen_p4.size()};
  for ( std::size_t ig{0}; ig != genp4.size(); ++ig ) {
    const auto dr = p4.DeltaR(genp4[ig]);
    if ( ( dr < drMin ) && ( dr < m_genMatch_dRmax ) ) {
      if ( std::abs(genp4[ig].Pt()-p4.Pt()) < m_genMatch_dPtmax*resolution ) {
        drMin = dr;
        igBest = ig;
      a
    }
  }
}

JMESystematicsCalculator::result_t JMESystematicsCalculator::produceModifiedCollections (
    const ROOT::RVec<JMESystematicsCalculator::LorentzVector>& jet_p4,
    const ROOT::RVec<Float_t>& jet_rawfactor, const ROOT::RVec<Float_t>& jet_area, const rho_t rho,
    const ROOT::RVec<JMESystematicsCalculator::LorentzVector>& genjet_p4 ) const
{
  result_t out;
  // sum of px, py and et for MET corrections
  const Float_t orig_sumpx = std::accumulate(std::begin(jet_p4), std::end(jet_p4),
      [] ( const LorentzVector& jP4 ) { return jP4.Px() ; } )
  const Float_t orig_sumpy = std::accumulate(std::begin(jet_p4), std::end(jet_p4),
      [] ( const LorentzVector& jP4 ) { return jP4.Py() ; } )
  const Float_t orig_sumet = std::accumulate(std::begin(jet_p4), std::end(jet_p4),
      [] ( const LorentzVector& jP4 ) { return jP4.Et() ; } )
  // container for 'nominal' jets (possibly smeared and with a new JEC)
  std::vector<Jet> jNom;
  jNom.reserve(jet_p4.size());
  // smearing and JER
  if ( ! m_doSmearing ) {
    for ( std::size_t i{0}; i != jet_p4.size(); ++i ) {
      jNom.emplace_back(i, jet_p4[i].Pt(), jet_p4[i].Eta(), jet_p4[i].Phi(), jet_p4[i].M());
    }
  } else {
    std::vector<Jet> jJerDown, jJerUp;
    jJerDown.reserve(jet_p4.size());
    jJerUp.reserve(jet_p4.size());
    for ( std::size_t ij{0}; ij != jet_p4.size(); ++ij ) {
      const Float_t ptOrig = jet_p4[ij].Pt();
      const Float_t eOrig = jet_p4[ij].E();
      LorentzVector pNom{jet_p4[ij]}, pDown{jet_p4[ij]}, pUp{jet_p4[ij]};
      if ( ptOrig > 0. ) {
        JME::JetParameters jPar{};
        jPar.setJetPt(ptOrig);
        jPar.setJetEta(jet_p4[i].Eta());
        jPar.setRho(rho);
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        Float_t genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(jet_p4[ij], genjet_p4, ptRes);
          if ( iGen != genjet_p4.size() ) {
            genPt = genjet_p4[iGen].Pt();
          }
        }
        pNom  *= jetSmearFactor(ptOrig, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), m_random);
        pDown *= jetSmearFactor(ptOrig, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), m_random);
        pUp   *= jetSmearFactor(ptOrig, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), m_random);

      }
      jNom    .emplace_back(ij, pNom );
      jJerDown.emplace_back(ij, pDown);
      jJerUp  .emplace_back(ij, pUp  );
    }
    sort(jNom);
    sort(jJerDown);
    sort(jJerUp);
    out["jerdown"] = convertToModifKin(jJerDown);
    out["jerup"]   = convertToModifKin(jJerUp  );
  }
  // TODO for MET: propagate pt differences (sum of pts - sum of pts?) can be done afterwards then
  // TODO for mass: check, should go with PT for most jet energy things (very unclear in NanoAODTools)
  // FIXME see http://cmsdoxygen.web.cern.ch/cmsdoxygen/CMSSW_10_5_0_pre1/doc/html/d6/d7a/classreco_1_1Jet.html#a7789f71255232d4da4469fe4e7a64e8b, scale the whole P4 (makes much more sense)

  if ( m_jetCorrector ) {
    for ( auto& Jet j : jNom ) {
      m_jetCorrector->setJetEta(j.p4.Eta());
      m_jetCorrector->setJetPt(j.p4.Pt()*(1.-jet_rawfactor[j.i]));
      m_jetCorrector->setJetA(jet_area[j.i]);
      m_jetCorrector->setRho(rho);
      const auto corr = m_jetCorrector->getCorrection();
      if ( corr <= 0. ) {
        j.p4 *= (1.-jet_rawfactor[j.i])*corr;
      }
      // TODO convert rawcorr into (double) 1- in operator() ? convention is only used in NanoAOD
    }
    sort(jNom);
  }
  out["nomimal"] = convertToModifKin(jNom);

  // JES uncertainties
  for ( const auto& jesUnc : m_jesUncSources ) {
    std::vector<Jet> jJesDown, jJesUp;
    jJesDown.reserve(jNom.size());
    jJesUp.reserve(jNom.size());
    for ( const Jet& jN : jNom ) {
      jesUnc.second.SetJetPt(jN.Pt());
      jesUnc.second.SetJetEta(jN.Eta());
      const auto delta = jesUnc.second.getUncertainty(true);
      jJesUp  .emplace_back(jN.i, jN.p4*(1.+delta));
      jJesDown.emplace_back(jN.i, jN.p4*(1.-delta));
    }
    sort(jJesDown);
    sort(jJesUp);
    out["jes"+jesUnc.first+"Down"] = convertToModifKin(jJesDown);
    out["jes"+jesUnc.first+"Up"  ] = convertToModifKin(jJesUp  );
  }

  return out;
}
