#include "JMESystematicsCalculators.h"

#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include "Math/VectorUtil.h"
#include "FactorizedJetCorrectorCalculator.h"

#include <cassert>

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

JetVariationsCalculator::result_t JetVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass )
{
  const auto nVariations = 1+( m_doSmearing ? 2 : 0 )+2*m_jesUncSources.size(); // 1(nom)+2(JER up/down)+2*len(JES)
  LogDebug << "JME:: hello from produce. Got " << jet_pt.size() << " jets" << std::endl;
  const auto nJets = jet_pt.size();
  result_t out{nVariations, jet_pt, jet_mass};
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
  std::size_t iVar = 1; // after nominal
  if ( m_doSmearing ) {
    LogDebug << "JME:: Smearing (seed=" << seed << ")" << std::endl;
    m_random.SetSeed(seed);
    p4compv_t pt_jerUp{pt_nom}, mass_jerUp{mass_nom};
    p4compv_t pt_jerDown{pt_nom}, mass_jerDown{mass_nom};
    for ( std::size_t i{0}; i != nJets; ++i ) {
      const float eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>(pt_nom[i], jet_phi[i], jet_eta[i], mass_nom[i]).E();
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
            LogDebug << "genPt=" << genPt << " ";
          }
        }
        const auto rand = ( genPt < 0. ) ? m_random.Gaus(0, ptRes) : -1.;
        LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
        const auto smearFactor_nom  = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
        const auto smearFactor_down = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
        const auto smearFactor_up   = jetSmearFactor(pt_nom[i], eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
        // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
        // LogDebug << "  smearfactors are NOMINAL=" << smearFactor_nom << ", DOWN=" << smearFactor_down << ", UP=" << smearFactor_up << std::endl;
        pt_nom[i]       *= smearFactor_nom;
        mass_nom[i]     *= smearFactor_nom;
        pt_jerDown[i]   *= smearFactor_down;
        mass_jerDown[i] *= smearFactor_down;
        pt_jerUp[i]     *= smearFactor_up;
        mass_jerUp[i]   *= smearFactor_up;
      }
    }
    out.set(iVar++, pt_jerUp  , mass_jerUp  );
    out.set(iVar++, pt_jerDown, mass_jerDown);
    LogDebug << "JME:: Done with smearing" << std::endl;
  } else {
    LogDebug << "JME:: No smearing" << std::endl;
  }
  out.set(0, pt_nom, mass_nom);

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
    out.set(iVar++, pt_jesUp  , mass_jesUp  );
    out.set(iVar++, pt_jesDown, mass_jesDown);
  }

#ifdef BAMBOO_JME_DEBUG
  assert(iVar == out.size());
  LogDebug << "JME:: returning " << out.size() << " modified jet collections" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: Jet_" << varNames[i] << ": ";
    for ( std::size_t j{0}; j != nJets; ++j ) {
      LogDebug << "(PT=" << out.pt(i)[j] << ", ETA=" << jet_eta[j] << ", PHI=" << jet_phi[j] << ", M=" << out.mass(i)[j] << ") ";
    }
    LogDebug << std::endl;
  }
#endif
  return out;
}

std::vector<std::string> JetVariationsCalculator::available() const
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


void Type1METVariationsCalculator::setL1JEC(const std::vector<JetCorrectorParameters>& jecParams)
{
  if ( ! jecParams.empty() ) {
    m_jetCorrectorL1 = std::unique_ptr<FactorizedJetCorrectorCalculator,jetcorrdeleter>{new FactorizedJetCorrectorCalculator(jecParams)};
  }
}

Type1METVariationsCalculator::result_t Type1METVariationsCalculator::produce(
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area,
    const p4compv_t& jet_muonSubtrFactor, const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF,
    const float rho,
    const std::uint32_t seed,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass,
    const float rawmet_phi, const float rawmet_pt,
    const float met_unclustenupdx, const float met_unclustenupdy,
    const p4compv_t& lowptjet_rawpt, const p4compv_t& lowptjet_eta, const p4compv_t& lowptjet_phi, const p4compv_t& lowptjet_area,
    const p4compv_t& lowptjet_muonSubtrFactor, const p4compv_t& lowptjet_neEmEF, const p4compv_t& lowptjet_chEmEF
    )
{
  const auto nVariations = 3+( m_doSmearing ? 3 : 0 )+2*m_jesUncSources.size(); // 1(nom)+2(unclust)+3(JER)+2*len(JES)
  result_t out{nVariations, rawmet_pt*std::cos(rawmet_phi), rawmet_pt*std::sin(rawmet_phi)};
  m_random.SetSeed(seed); // just in case
  LogDebug << "JME:: hello from produce. Got " << jet_pt.size() << " jets and " << lowptjet_rawpt.size() << " low-PT jets" << std::endl;
  // normal jets
  addVariations(out, jet_pt, jet_eta, jet_phi, jet_mass,
      jet_rawcorr, jet_area, jet_muonSubtrFactor, jet_neEmEF, jet_chEmEF,
      ROOT::VecOps::RVec<bool>(), rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass);
  // low-PT jets
  p4compv_t lowptjet_zero(lowptjet_rawpt.size(), 0.);
  addVariations(out, lowptjet_rawpt, lowptjet_eta, lowptjet_phi, lowptjet_zero,
      lowptjet_zero, lowptjet_area, lowptjet_muonSubtrFactor,
      ( lowptjet_neEmEF.empty() ? lowptjet_zero : lowptjet_neEmEF  ),
      ( lowptjet_chEmEF.empty() ? lowptjet_zero : lowptjet_chEmEF  ),
      ROOT::VecOps::RVec<bool>(), rho, genjet_pt, genjet_eta, genjet_phi, genjet_mass);
  // unclustered energy, based on nominal (0)
  out.setXY(nVariations-2, out.px(0)+met_unclustenupdx, out.py(0)+met_unclustenupdy);
  out.setXY(nVariations-1, out.px(0)-met_unclustenupdx, out.py(0)-met_unclustenupdy);

#ifdef BAMBOO_JME_DEBUG
  LogDebug << "JME:: returning " << out.size() << " modified METs" << std::endl;
  const auto varNames = available();
  assert(varNames.size() == nVariations);
  for ( std::size_t i{0}; i != nVariations; ++i ) {
    LogDebug << "JME:: MET_" << varNames[i] << ": PT=" << out.pt(i) << ", PHI=" << out.phi(i) << std::endl;
  }
#endif
  return out;
}

// for a single jet collection
void Type1METVariationsCalculator::addVariations(Type1METVariationsCalculator::result_t& out,
    const p4compv_t& jet_pt, const p4compv_t& jet_eta, const p4compv_t& jet_phi, const p4compv_t& jet_mass,
    const p4compv_t& jet_rawcorr, const p4compv_t& jet_area, const p4compv_t& jet_muonSubtrFactor,
    const p4compv_t& jet_neEmEF, const p4compv_t& jet_chEmEF, const ROOT::VecOps::RVec<bool>& jet_mask, const float rho,
    const p4compv_t& genjet_pt, const p4compv_t& genjet_eta, const p4compv_t& genjet_phi, const p4compv_t& genjet_mass
    )
{
  FactorizedJetCorrectorCalculator::VariableValues vals, valsL1;
  const auto nJets = jet_pt.size();
  for ( std::size_t i{0}; i != nJets; ++i ) {
    if ( jet_mask.empty() || jet_mask[i] ) {
      // L1 and full (L1L2L3) JEC for muon-subtracted jet
      auto jet_pt_nom = jet_pt[i];
      auto jet_mass_nom = jet_mass[i];
      const auto jet_pt_raw = jet_pt[i]*(1-jet_rawcorr[i]);
      vals.setJetEta(jet_eta[i]);
      vals.setJetPt(jet_pt_raw);
      vals.setJetA(jet_area[i]);
      vals.setRho(rho);
      LogDebug << "Jet #" << i << " ETA=" << jet_eta[i] << ", PT_raw=" << jet_pt_raw << ", area=" << jet_area[i] << std::endl;
      auto jecL1L2L3 = m_jetCorrector->getCorrection(vals);
      if ( jecL1L2L3 <= 0. ) {
        jecL1L2L3 = 1.;
      } else {
        jet_pt_nom = jet_pt_raw*jecL1L2L3;
        jet_mass_nom = jet_mass_nom*(1-jet_rawcorr[i])*jecL1L2L3;
      }
      valsL1.setJetEta(jet_eta[i]);
      valsL1.setJetPt(jet_pt_raw);
      valsL1.setJetA(jet_area[i]);
      valsL1.setRho(rho);
      auto jecL1 = m_jetCorrectorL1->getCorrection(valsL1);
      if ( jecL1     <= 0. ) { jecL1     = 1.; }
      const float jet_pt_raw_nomu = jet_pt_raw*(1-jet_muonSubtrFactor[i]);
      const float muon_pt = jet_pt_raw*jet_muonSubtrFactor[i];
      float jet_pt_nomuL1L2L3{jet_pt_raw_nomu}, jet_pt_nomuL1{jet_pt_raw_nomu};
      if ( jet_pt_raw_nomu*jecL1L2L3 > m_unclEnThreshold ) {
        jet_pt_nomuL1L2L3 = jet_pt_raw_nomu*jecL1L2L3;
        jet_pt_nomuL1     = jet_pt_raw_nomu*jecL1;
      }
      const auto jet_pt_L1L2L3 = jet_pt_nomuL1L2L3 + muon_pt;
      const auto jet_pt_L1     = jet_pt_nomuL1     + muon_pt;
      LogDebug << "jecL1L2L3=" << jecL1L2L3 << ", jecL1=" << jecL1 << "; PT_L1L2L3=" << jet_pt_L1L2L3 << ", PT_L1=" << jet_pt_L1 << ", PT_mu=" << muon_pt << std::endl;
      // JER / smearing
      float jerF_nom{1.}, jerF_up{1.}, jerF_down{1.};
      if ( m_doSmearing ) {
        const float eOrig = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>(jet_pt_nom, jet_phi[i], jet_eta[i], jet_mass_nom).E();
        if ( jet_pt_nom > 0. ) {
          JME::JetParameters jPar{
              {JME::Binning::JetPt , jet_pt_nom},
              {JME::Binning::JetEta, jet_eta[i]},
              {JME::Binning::Rho   , rho} };
          const auto ptRes  = m_jetPtRes.getResolution(jPar);
          LogDebug << "JME:: JetParameters: pt=" << jet_pt_nom << ", eta=" << jet_eta[i] << ", rho=" << rho << "; ptRes=" << ptRes << std::endl;
          LogDebug << "JME:: ";
          float genPt = -1;
          if ( m_smearDoGenMatch ) {
            const auto iGen = findGenMatch(jet_pt_nom, jet_eta[i], jet_phi[i], genjet_pt, genjet_eta, genjet_phi, ptRes*jet_pt_nom);
            if ( iGen != genjet_pt.size() ) {
              genPt = genjet_pt[iGen];
              LogDebug << "genPt=" << genPt << " ";
            }
          }
          const auto rand = ( genPt < 0. ) ? m_random.Gaus(0, ptRes) : -1.;
          LogDebug << "jet_pt_resolution: " << ptRes << ", rand: " << rand << std::endl;
          jerF_nom  = jetSmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL), rand);
          jerF_down = jetSmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::DOWN   ), rand);
          jerF_up   = jetSmearFactor(jet_pt_nom, eOrig, genPt, ptRes, m_jetEResSF.getScaleFactor(jPar, Variation::UP     ), rand);
          // LogDebug << "  scalefactors are NOMINAL=" << m_jetEResSF.getScaleFactor(jPar, Variation::NOMINAL) << ", DOWN=" << m_jetEResSF.getScaleFactor(jPar, Variation::DOWN) << ", UP=" << m_jetEResSF.getScaleFactor(jPar, Variation::UP) << std::endl;
          // LogDebug << "  smearfactors are NOMINAL=" << jerF_nom << ", DOWN=" << jerF_down << ", UP=" << jerF_up << std::endl;
        }
      }
      if ( ( jet_pt_L1L2L3 > m_unclEnThreshold ) && ( (jet_neEmEF[i]+jet_chEmEF[i]) < 0.9 ) ) {
        std::size_t iVar = 0;
        const auto jet_cosPhi = std::cos(jet_phi[i]);
        const auto jet_sinPhi = std::sin(jet_phi[i]);
        out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3);             // nominal
        if ( m_doSmearing ) {
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3*jerF_nom);  // JER
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3*jerF_up);   // JER-up
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3*jerF_down); // JER-down
        }
        for ( auto& jesUnc : m_jesUncSources ) {
          jesUnc.second.setJetPt(jet_pt_L1L2L3);
          jesUnc.second.setJetEta(jet_eta[i]);
          const auto delta = jesUnc.second.getUncertainty(true);
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3*(1+delta)); // JES_i-up
          out.addR_proj(iVar++, jet_cosPhi, jet_sinPhi, jet_pt_L1 - jet_pt_L1L2L3*(1-delta)); // JES_i-down
        }
#ifdef BAMBOO_JME_DEBUG
        assert(iVar+2 == out.size()); // last two are unclustered energy up and down
#endif
      }
    }
  }
}

std::vector<std::string> Type1METVariationsCalculator::available() const
{
  std::vector<std::string> products = { "nominal" };
  if ( m_doSmearing ) {
    products.emplace_back("jer");
    products.emplace_back("jerup");
    products.emplace_back("jerdown");
  }
  for ( const auto& src : m_jesUncSources ) {
    products.emplace_back("jes"+src.first+"up");
    products.emplace_back("jes"+src.first+"down");
  }
  products.emplace_back("unclustEnup");
  products.emplace_back("unclustEndown");
  return products;
}
