#include "RochesterCorrectionCalculator.h"

#include "RoccoR.h"

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

RochesterCorrectionCalculator::result_t RochesterCorrectionCalculator::produceModifiedCollections(
      const p4compv_t& muon_pt, const p4compv_t& muon_eta, const p4compv_t& muon_phi, const p4compv_t& muon_mass, const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const ROOT::VecOps::RVec<Int_t>& muon_genIdx, const p4compv_t& gen_pt) const
{
  std::vector<Muon> muons;
  muons.reserve(muon_pt.size());
  for ( std::size_t im{0}; im != muon_pt.size(); ++im ) {
    muons.emplace_back(im, muon_pt[im], muon_eta[im], muon_phi[im], muon_mass[im]);
  }
  p4compv_t muonGenpt;
  if ( m_roccor && ( ! muon_genIdx.empty() ) && ( ! gen_pt.empty() ) ) {
    for ( std::size_t i{0}; i < muon_genIdx.size(); ++i ) {
      muonGenpt.push_back(muon_genIdx[i] != -1 ? gen_pt[muon_genIdx[i]] : -1.);
    }
  }
  return produceModifiedCollections(muons, muon_charge, muon_nlayers, muonGenpt);
}

namespace {

RochesterCorrectionCalculator::result_entry_t convertToModifKin( const std::vector<RochesterCorrectionCalculator::Muon>& muons )
{
  RochesterCorrectionCalculator::result_entry_t::Indices idx;
  RochesterCorrectionCalculator::result_entry_t::Momenta mom;
  idx.reserve(muons.size());
  mom.reserve(muons.size());
  for ( const auto& j : muons ) {
    // TODO copies can be made implicit when push_back(const&) is there
    idx.push_back(std::size_t{j.i});
    mom.push_back(RochesterCorrectionCalculator::LorentzVector{j.p4});
  }
  return RochesterCorrectionCalculator::result_entry_t(std::move(idx), std::move(mom));
}

void sort( std::vector<RochesterCorrectionCalculator::Muon>& muons )
{
  std::sort(std::begin(muons), std::end(muons), [] ( const RochesterCorrectionCalculator::Muon& m1, const RochesterCorrectionCalculator::Muon& m2 ) { return m1.p4.Pt() > m2.p4.Pt(); });
}

}

RochesterCorrectionCalculator::result_t RochesterCorrectionCalculator::produceModifiedCollections(std::vector<RochesterCorrectionCalculator::Muon>& muons, const ROOT::VecOps::RVec<Int_t>& muon_charge, const ROOT::VecOps::RVec<Int_t>& muon_nlayers, const p4compv_t& muon_genpt) const
{
  LogDebug << "Rochester:: hello from produceModifiedCollections. Got " << muons.size() << " muons" << std::endl;
  result_t out;
  if ( ! m_roccor ) {
    LogDebug << "Rochester:: No correction" << std::endl;
    out["nominal"] = convertToModifKin(muons);
  } else {
    if ( muon_genpt.empty() ) { // DATA
      for ( std::size_t i{0}; i < muons.size(); ++i ) {
        const auto sf = m_roccor->kScaleDT(muon_charge[i], muons[i].p4.Pt(), muons[i].p4.Eta(), muons[i].p4.Phi());
        LogDebug << "Rochester:: DATA scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muons[i].p4.Pt() << " ETA=" << muons[i].p4.Eta() << " PHI=" << muons[i].p4.Phi() << std::endl;
        muons[i].p4 *= sf;
      }
    } else { // MC
      for ( std::size_t i{0}; i < muons.size(); ++i ) {
        if ( muon_genpt[i] > -1. ) { // gen muon available
          const auto sf = m_roccor->kSpreadMC(muon_charge[i], muons[i].p4.Pt(), muons[i].p4.Eta(), muons[i].p4.Phi(), muon_genpt[i]);
          LogDebug << "Rochester:: MC scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muons[i].p4.Pt() << " ETA=" << muons[i].p4.Eta() << " PHI=" << muons[i].p4.Phi() << " GENPT=" << muon_genpt[i] << " (nLayers=" << muon_nlayers[i] << ")" << std::endl;
          muons[i].p4 *= sf;
        } else { // gen muon not available
          std::uniform_real_distribution<> d(0., 1.);
          const auto sf = m_roccor->kSmearMC(muon_charge[i], muons[i].p4.Pt(), muons[i].p4.Eta(), muons[i].p4.Phi(), muon_nlayers[i], d(m_random));
          LogDebug << "Rochester:: MC scale " << sf << " for charge=" << muon_charge[i] << " PT=" << muons[i].p4.Pt() << " ETA=" << muons[i].p4.Eta() << " PHI=" << muons[i].p4.Phi() << " NLayers=" << muon_nlayers[i] << std::endl;
          muons[i].p4 *= sf;
        }
      }
    }
    sort(muons);
    out["nominal"] = convertToModifKin(muons);
  }
  return out;
}

bool RochesterCorrectionCalculator::hasProduct(const std::string& name) const
{
  if ( name == "nominal" )
    return true;
  return false;
}