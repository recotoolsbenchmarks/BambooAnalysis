#include "JMESystematicsCalculator.h"

#include "Math/VectorUtil.h"
#include "FactorizedJetCorrectorCalculator.h"

// #define BAMBOO_JME_DEBUG // uncomment to debug

#ifdef BAMBOO_JME_DEBUG
#define LogDebug std::cout
#else
#define LogDebug if (false) std::cout
#endif

JetMETVariationsCalculatorBase::~JetMETVariationsCalculatorBase()
{}

void JetMETVariationsCalculatorBase::setJEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrector = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

void JetMETVariationsCalculatorBase::jetcorrdeleter::operator() (FactorizedJetCorrectorCalculator* ptr) const
{ delete ptr; }

namespace {
  float jetSmearFactor( float pt, float eOrig, float genPt, float ptRes, float sfUncert, float rand )
  {
    float smear = 1.;
    if ( genPt > 0. ) {
      smear = 1. + (sfUncert-1.)*(pt - genPt)/pt;
    } else if ( sfUncert > 1. ) {
      smear = 1. + rand*std::sqrt(sfUncert*sfUncert-1.);
    }
    if ( smear*eOrig < 1.e-2 ) {
      smear = 1.e-2/eOrig;
    }
    return smear;
  }

  // because something goes wrong with linking ROOT::Math::VectorUtil::Phi_mpi_pi
  template<typename T>
  T phi_mpi_pi(T angle) {
    if ( angle <= M_PI && angle > -M_PI ) {
      return angle;
    }
    if ( angle > 0 ) {
      const int n = static_cast<int>(.5*(angle*M_1_PI+1.));
      angle -= 2*n*M_PI;
    } else {
      const int n = static_cast<int>(-.5*(angle*M_1_PI-1.));
      angle += 2*n*M_PI;
    }
    return angle;
  }
}

// TODO with orig MET and jets (sumpx,sumpy): calc modif MET(sig), produce bigger results type

std::size_t JetMETVariationsCalculatorBase::findGenMatch(const float pt, const float eta, const float phi, const ROOT::VecOps::RVec<float>& gen_pt, const ROOT::VecOps::RVec<float>& gen_eta, const ROOT::VecOps::RVec<float>& gen_phi, const float resolution ) const
{
  auto dr2Min = std::numeric_limits<float>::max();
  std::size_t igBest{gen_pt.size()};
  LogDebug << "(DRs: ";
  for ( std::size_t ig{0}; ig != gen_pt.size(); ++ig ) {
    const auto dphi = phi_mpi_pi(gen_phi[ig]-phi);
    const auto deta = (gen_eta[ig]-eta);
    const auto dr2 = dphi*dphi + deta*deta;
    LogDebug << "dr2=" << dr2;
    if ( ( dr2 < dr2Min ) && ( dr2 < m_genMatch_dR2max ) ) {
      LogDebug << "->dpt=" << std::abs(gen_pt[ig]-pt) << ",res=" << resolution;
      if ( std::abs(gen_pt[ig]-pt) < m_genMatch_dPtmax*resolution ) {
        LogDebug << "->best:" << ig;
        dr2Min = dr2;
        igBest = ig;
      }
    }
    LogDebug << ", ";
  }
  LogDebug << ")";
  return igBest;
}

JMESystematicsCalculator::result_t JMESystematicsCalculator::produceModifiedCollections(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass )
{
  LogDebug << "JME:: hello from produceModifiedCollections. Got " << jet_pt.size() << " jets" << std::endl;
  const auto nJets = jet_pt.size();
  result_t out{jet_eta, jet_phi};
  p4compv_t pt_nom{jet_pt}, mass_nom{jet_mass};
  if ( m_jetCorrector ) {
    LogDebug << "JME:: reapplying JEC" << std::endl;
    FactorizedJetCorrectorCalculator::VariableValues vals;
    for ( std::size_t i{0}; i != nJets; ++i ) {
      vals.setJetEta(jet_eta[i]);
      vals.setJetPt(jet_pt[i]*(1.-jet_rawcorr[i]));
      vals.setJetA(jet_area[i]);
      vals.setRho(rho);
      const auto corr = m_jetCorrector->getCorrection(vals);
      if ( corr > 0. ) {
        const auto newc = (1.-jet_rawcorr[i])*corr;
        pt_nom[i]   *= newc;
        mass_nom[i] *= newc;
      }
    }
#ifdef BAMBOO_JME_DEBUG
    LogDebug << "JME:: with reapplied JEC: ";
    for ( std::size_t i{0}; i != nJets; ++i ) {
      LogDebug << "(PT=" << pt_nom[i] << ", ETA=" << jet_eta[i] << ", PHI=" << jet_phi[i] << ", M=" << mass_nom[i] << ") ";
    }
    LogDebug << std::endl;
#endif
  } else {
    LogDebug << "JME:: Not reapplying JEC" << std::endl;
  }
  // smearing and JER
  if ( m_doSmearing ) {
    LogDebug << "JME:: Smearing (seed=" << seed << ")" << std::endl;
    m_random.SetSeed(seed);
    p4compv_t pt_jerUp{pt_nom}, mass_jerUp{mass_nom};
    p4compv_t pt_jerDown{pt_nom}, mass_jerDown{mass_nom};
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const float eOrig = LorentzVector(pt_nom[i], jet_phi[i], jet_eta[i], mass_nom[i]).E();
      if ( pt_nom[i] > 0. ) {
        JME::JetParameters jPar{
            {JME::Binning::JetPt , pt_nom[i]},
            {JME::Binning::JetEta, jet_eta[i]},
            {JME::Binning::Rho   , rho} };
        const auto ptRes  = m_jetPtRes.getResolution(jPar);
        LogDebug << "JME:: JetParameters: pt=" << pt_nom[i] << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
        LogDebug << "JME:: ";
        float genPt = -1;
        if ( m_smearDoGenMatch ) {
          const auto iGen = findGenMatch(pt_nom[i], jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*pt_nom[i]);
          if ( iGen != genjet_pt.size() ) {
            genPt = genjet_pt[iGen];
            LogDebug << "genPt=" << genPt;
          }
        }
        const auto rand = ( genPt < 0. ) ? m_random.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        const auto smearFactor_nom  = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        const auto smearFactor_down = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        const auto smearFactor_up   = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        LogDebug << "  smearfactors are NOMINAL=" << smearFactor_nom << ", DOWN=" << smearFactor_down << ", UP=" << smearFactor_up << std::endl;
        pt_nom[i]       *= smearFactor_nom;
        mass_nom[i]     *= smearFactor_nom;
        pt_jerDown[i]   *= smearFactor_down;
        mass_jerDown[i] *= smearFactor_down;
        pt_jerUp[i]     *= smearFactor_up;
        mass_jerUp[i]   *= smearFactor_up;
      }
    }
    out.add("jerdown", pt_jerDown, mass_jerDown);
    out.add("jerup"  , pt_jerUp  , mass_jerUp  );
    LogDebug << "JME:: Done with smearing" << std::endl;
  } else {
    LogDebug << "JME:: No smearing" << std::endl;
  }
  out.add("nominal", pt_nom, mass_nom);

  // JES uncertainties
  for ( auto& jesUnc : m_jesUncSources ) {
    LogDebug << "JME:: evaluating JES uncertainty: " << jesUnc.first << std::endl;
    p4compv_t pt_jesDown{pt_nom}, mass_jesDown{mass_nom};
    p4compv_t pt_jesUp{pt_nom}, mass_jesUp{mass_nom};
    for ( std::size_t i{0}; i != nJets; ++i ) {
      jesUnc.second.setJetPt(pt_nom[i]);
      jesUnc.second.setJetEta(jet_eta[i]);
      const auto delta = jesUnc.second.getUncertainty(true);
      pt_jesDown[i]   *= (1.-delta);
      mass_jesDown[i] *= (1.-delta);
      pt_jesUp[i]     *= (1.+delta);
      mass_jesUp[i]   *= (1.+delta);
    }
    out.add("jes"+jesUnc.first+"down", pt_jesDown, mass_jesDown);
    out.add("jes"+jesUnc.first+"up"  , pt_jesUp  , mass_jesUp  );
  }

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  for ( const auto& entryName : availableProducts() ) {
    LogDebug << "JME:: " << entryName << ": ";
    const auto& entry = out[entryName];
    for ( std::size_t i{0}; i != nJets; ++i ) {
      LogDebug << "(PT=" << entry.pt()[i] << ", ETA=" << entry.eta()[i] << ", PHI=" << entry.phi()[i] << ", M=" << entry.mass()[i] << ") ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

std::vector<std::string> JMESystematicsCalculator::availableProducts() const
{
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    products.emplace_back("jerup");
    products.emplace_back("jerdown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  return products;
}
