#include "JMESystematicsCalculator.h"

#include "Math/VectorUtil.h"
#include "FactorizedJetCorrectorCalculator.h"

#ifdef BAMBOO_JME_DEBUG
#define LogDebug std::cout
#else
#define LogDebug if (false) std::cout
#endif

JMESystematicsCalculator::~JMESystematicsCalculator()
{}

void JMESystematicsCalculator::setJEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrector = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

void JMESystematicsCalculator::jetcorrdeleter::operator() (FactorizedJetCorrectorCalculator* ptr) const
{ delete ptr; }

JMESystematicsCalculator::result_t JMESystematicsCalculator::produceModifiedCollections(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const float rho,
    const float met_phi, const float met_pt, const float sumEt,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass )
{
  ROOT::VecOps::RVec<LorentzVector> jetP4;
  jetP4.reserve(jet_pt.size());
  for ( std::size_t i{0}; i < jet_pt.size(); ++i ) {
    jetP4.emplace_back(jet_pt[i], jet_eta[i], jet_phi[i], jet_mass[i]);
  }
  ROOT::VecOps::RVec<LorentzVector> genP4;
  genP4.reserve(genjet_pt.size());
  for ( std::size_t i{0}; i < genjet_pt.size(); ++i ) {
    genP4.emplace_back(genjet_pt[i], genjet_eta[i], genjet_phi[i], genjet_mass[i]);
  }
  ROOT::VecOps::RVec<double> jet_rawfactor{};
  jet_rawfactor.reserve(jet_rawcorr.size());
  for ( const auto rawcorr : jet_rawcorr ) {
    jet_rawfactor.push_back(1.-rawcorr);
  }
  return produceModifiedCollections(jetP4, jet_rawfactor, jet_area, rho, genP4);
}

namespace {
  float jetSmearFactor( float pt, float eOrig, float genPt, float ptRes, float sfUncert, std::mt19937& randgen )
  {
    float smear = 1.;
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
}

// TODO with orig MET and jets (sumpx,sumpy): calc modif MET(sig), produce bigger results type
// TODO non-const ref and include sort then
JMESystematicsCalculator::result_entry_t JMESystematicsCalculator::convertToModifKin( const std::vector<JMESystematicsCalculator::Jet>& jets ) const
{
  result_entry_t::Indices idx;
  result_entry_t::Momenta mom;
  idx.reserve(jets.size());
  mom.reserve(jets.size());
  for ( const auto& j : jets ) {
    // TODO copies can be made implicit when push_back(const&) is there
    idx.push_back(std::size_t{j.i});
    mom.push_back(LorentzVector{j.p4});
  }
  return JMESystematicsCalculator::result_entry_t(std::move(idx), std::move(mom));
}

void JMESystematicsCalculator::sort( std::vector<Jet>& jets ) const
{
  std::sort(std::begin(jets), std::end(jets), [] ( const Jet& j1, const Jet& j2 ) { return j1.p4.Pt() > j2.p4.Pt(); });
}

std::size_t JMESystematicsCalculator::findGenMatch( const JMESystematicsCalculator::LorentzVector& p4, const ROOT::VecOps::RVec<JMESystematicsCalculator::LorentzVector>& genp4, const float resolution ) const
{
  auto drMin = std::numeric_limits<float>::max();
  std::size_t igBest{genp4.size()};
  LogDebug << "(DRs: ";
  for ( std::size_t ig{0}; ig != genp4.size(); ++ig ) {
    const auto dr = ROOT::Math::VectorUtil::DeltaR(p4, genp4[ig]);
    LogDebug << dr;
    if ( ( dr < drMin ) && ( dr < m_genMatch_dRmax ) ) {
      LogDebug << "->dpt=" << std::abs(genp4[ig].Pt()-p4.Pt()) << ",res=" << resolution;
      if ( std::abs(genp4[ig].Pt()-p4.Pt()) < m_genMatch_dPtmax*resolution ) {
        LogDebug << "->best:" << ig;
        drMin = dr;
        igBest = ig;
      }
    }
    LogDebug << ", ";
  }
  LogDebug << ")";
  return igBest;
}

JMESystematicsCalculator::result_t JMESystematicsCalculator::produceModifiedCollections (
    const ROOT::VecOps::RVec<JMESystematicsCalculator::LorentzVector>& jet_p4,
    const ROOT::VecOps::RVec<double>& jet_rawfactor, const JMESystematicsCalculator::p4compv_t& jet_area, float rho,
    const ROOT::VecOps::RVec<JMESystematicsCalculator::LorentzVector>& genjet_p4 )
{
  LogDebug << "JME:: hello from produceModifiedCollections. Got " << jet_p4.size() << " jets" << std::endl;
  result_t out;
  // sum of px, py and et for MET corrections
  const float orig_sumpx = std::accumulate(std::begin(jet_p4), std::end(jet_p4), 0.,
      [] ( float res, const LorentzVector& jP4 ) { return res+jP4.Px() ; } );
  const float orig_sumpy = std::accumulate(std::begin(jet_p4), std::end(jet_p4), 0.,
      [] ( float res, const LorentzVector& jP4 ) { return res+jP4.Py() ; } );
  const float orig_sumet = std::accumulate(std::begin(jet_p4), std::end(jet_p4), 0.,
      [] ( float res, const LorentzVector& jP4 ) { return res+jP4.Et() ; } );
  // container for 'nominal' jets (possibly smeared and with a new JEC)
  LogDebug << "JME:: MET sum (px,py)=(" << orig_sumpx << ", " << orig_sumpy << ") sumET=" << orig_sumet << std::endl;
  std::vector<Jet> jNom;
  jNom.reserve(jet_p4.size());
  // smearing and JER
  if ( ! m_doSmearing ) {
    for ( std::size_t i{0}; i != jet_p4.size(); ++i ) {
      jNom.emplace_back(i, jet_p4[i].Pt(), jet_p4[i].Eta(), jet_p4[i].Phi(), jet_p4[i].M());
    }
    LogDebug << "JME:: No smearing" << std::endl;
  } else {
    LogDebug << "JME:: Smearing" << std::endl;
    std::vector<Jet> jJerDown, jJerUp;
    jJerDown.reserve(jet_p4.size());
    jJerUp.reserve(jet_p4.size());
    for ( std::size_t ij{0}; ij != jet_p4.size(); ++ij ) {
      const float ptOrig = jet_p4[ij].Pt();
      const float eOrig = jet_p4[ij].E();
      LorentzVector pNom{jet_p4[ij]}, pDown{jet_p4[ij]}, pUp{jet_p4[ij]};
      if ( ptOrig > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , ptOrig},
            {JME::Binning::JetEta, jet_p4[ij].Eta()},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(jet_p4[ij], genjet_p4, ptRes*ptOrig);
          if ( iGen != genjet_p4.size() ) {
            genPt = genjet_p4[iGen].Pt();
            LogDebug << "genPt=" << genPt;
          }
        }
        LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
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
    LogDebug << "JME:: Done with smearing" << std::endl;
  }
  // TODO for MET: propagate pt differences (sum of pts - sum of pts?) can be done afterwards then
  if ( m_jetCorrector ) {
    LogDebug << "JME:: reapplying JEC" << std::endl;
    FactorizedJetCorrectorCalculator::VariableValues vals;
    for ( auto& j : jNom ) {
      vals.setJetEta(j.p4.Eta());
      vals.setJetPt(j.p4.Pt()*jet_rawfactor[j.i]);
      vals.setJetA(jet_area[j.i]);
      vals.setRho(rho);
      const auto corr = m_jetCorrector->getCorrection(vals);
      if ( corr <= 0. ) {
        j.p4 *= jet_rawfactor[j.i]*corr;
      }
    }
    sort(jNom);
  }
  out["nomimal"] = convertToModifKin(jNom);

  // JES uncertainties
  for ( auto& jesUnc : m_jesUncSources ) {
    LogDebug << "JME:: evaluating JES uncertainty: " << jesUnc.first << std::endl;
    std::vector<Jet> jJesDown, jJesUp;
    jJesDown.reserve(jNom.size());
    jJesUp.reserve(jNom.size());
    for ( const Jet& jN : jNom ) {
      jesUnc.second.setJetPt(jN.p4.Pt());
      jesUnc.second.setJetEta(jN.p4.Eta());
      const auto delta = jesUnc.second.getUncertainty(true);
      jJesUp  .emplace_back(jN.i, jN.p4*(1.+delta));
      jJesDown.emplace_back(jN.i, jN.p4*(1.-delta));
    }
    sort(jJesDown);
    sort(jJesUp);
    out["jes"+jesUnc.first+"Down"] = convertToModifKin(jJesDown);
    out["jes"+jesUnc.first+"Up"  ] = convertToModifKin(jJesUp  );
  }

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  for ( const auto& entry : out ) {
    LogDebug << "JME:: " << entry.first << ": ";
    for ( std::size_t i{0}; i != entry.second.momenta().size(); ++i ) {
      const auto mom = entry.second.momenta()[i];
      LogDebug << "(PT=" << mom.Pt() << ", ETA=" << mom.Eta() << ", PHI=" << mom.Phi() << ", M=" << mom.M() << ")" << entry.second.indices()[i] << "  ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

bool JMESystematicsCalculator::hasProduct(const std::string& name) const
{
  if ( name == "nominal" )
    return true;
  else if ( m_doSmearing && ( ( name == "jerup" ) || ( name == "jerdown" ) ) )
    return true;
  else if ( ( ! m_jesUncSources.empty() ) && ( name.substr(0,3) == "jes" ) ) {
    std::string jesunc{};
    if ( name.substr(name.size()-4) == "Down" ) {
      jesunc = name.substr(3, name.size()-4);
    } else if ( name.substr(name.size()-2) == "Up" ) {
      jesunc = name.substr(3, name.size()-2);
    }
    return ( ! jesunc.empty() ) && ( std::end(m_jesUncSources) != m_jesUncSources.find(jesunc) );
  }
  return false;
}
