#include "RochesterCorrectionCalculator.h"

#include "RoccoR.h"

#include "bamboorandom.h"

// #define BAMBOO_ROCCOR_DEBUG // uncomment to debug

#ifdef BAMBOO_ROCCOR_DEBUG
#define LogDebug std::cout
#else
#define LogDebug if (false) std::cout
#endif

RochesterCorrectionCalculator::~RochesterCorrectionCalculator()
{}

void RochesterCorrectionCalculator::setRochesterCorrection(const std::string& params)
{
  if ( ! params.empty() ) {
    m_roccor = std::unique_ptr<RoccoR,roccordeleter>{new RoccoR(params)};
  }
}

void RochesterCorrectionCalculator::roccordeleter::operator() (RoccoR* ptr) const
{ delete ptr; }

RochesterCorrectionCalculator::result_t RochesterCorrectionCalculator::produce(
      const p4compv_t& muon_pt, const p4compv_t& muon_eta, const p4compv_t& muon_phi, const p4compv_t& muon_mass, const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const ROOT::VecOps::RVec<Int_t>& muon_genIdx, const p4compv_t& gen_pt, const uint32_t seed) const
{
  const auto nVariations = 1;
  LogDebug << "Rochester:: hello from produce. Got " << muon_pt.size() << " muons" << std::endl;
  result_t out{nVariations, muon_pt};
  if ( ! m_roccor ) {
    LogDebug << "Rochester:: No correction" << std::endl;
  } else {
    p4compv_t nom_pt = muon_pt;
    if ( muon_genIdx.empty() ) { // DATA
      for ( std::size_t i{0}; i < muon_pt.size(); ++i ) {
        const auto sf = m_roccor->kScaleDT(muon_charge[i], muon_pt[i], muon_eta[i], muon_phi[i]);
        LogDebug << "Rochester:: DATA scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muon_pt[i] << " ETA=" << muon_eta[i] << " PHI=" << muon_phi[i] << std::endl;
        nom_pt[i] = muon_pt[i]*sf;
      }
    } else { // MC
      auto& rg = rdfhelpers::getStdMT19937(seed);
      for ( std::size_t i{0}; i < muon_pt.size(); ++i ) {
        if ( muon_genIdx[i] != -1 ) { // gen muon available
          const auto sf = m_roccor->kSpreadMC(muon_charge[i], muon_pt[i], muon_eta[i], muon_phi[i], gen_pt[muon_genIdx[i]]);
          LogDebug << "Rochester:: MC scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muon_pt[i] << " ETA=" << muon_eta[i] << " PHI=" << muon_phi[i] << " GENPT=" << gen_pt[muon_genIdx[i]] << " (nLayers=" << muon_nlayers[i] << ")" << std::endl;
          nom_pt[i] = muon_pt[i]*sf;
        } else { // gen muon not available
          std::uniform_real_distribution<> d(0., 1.);
          const auto sf = m_roccor->kSmearMC(muon_charge[i], muon_pt[i], muon_eta[i], muon_phi[i], muon_nlayers[i], d(rg));
          LogDebug << "Rochester:: MC scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muon_pt[i] << " ETA=" << muon_eta[i] << " PHI=" << muon_phi[i] << " NLayers=" << muon_nlayers[i] << std::endl;
          nom_pt[i] = muon_pt[i]*sf;
        }
      }
    }
    out.set(0, nom_pt);
  }
  return out;
}

std::vector<std::string> RochesterCorrectionCalculator::available() const
{
  return { "nominal" }; // TODO add systematics
}
