#pragma once

#include <ROOT/RVec.hxx>
#include <Math/GenVector/LorentzVector.h>
#include <Math/GenVector/PtEtaPhiM4D.h>
#include <boost/container/flat_map.hpp>

namespace rdfhelpers {

class ModifiedPtCollection { // map of variation collections
public:
  using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
  using p4compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtCollection() = default;
  ModifiedPtCollection(p4compv_t eta, p4compv_t phi, p4compv_t mass) : m_eta(eta), m_phi(phi), m_mass(mass) {}
  void add(const std::string& name, p4compv_t pt)
  {
    m_variations[name] = pt;
  }

  const p4compv_t& eta() const { return m_eta; }
  const p4compv_t& phi() const { return m_phi; }
  const p4compv_t& mass() const { return m_mass; }

  class Variation { // single variation
  public:
    Variation(const ModifiedPtCollection& parent, const p4compv_t& pt) : m_parent(parent), m_pt(pt) {}
    const p4compv_t& pt  () const { return m_pt; }
    const p4compv_t& eta () const { return m_parent.eta(); }
    const p4compv_t& phi () const { return m_parent.phi(); }
    const p4compv_t& mass() const { return m_parent.mass(); }
    std::size_t size() const { return m_pt.size(); }
  private:
    const ModifiedPtCollection& m_parent;
    const p4compv_t& m_pt;
  };

  Variation at(const std::string& name) const { return Variation(*this, m_variations.at(name)); }
  std::size_t size() const { return m_variations.size(); }
private:
  boost::container::flat_map<std::string,p4compv_t> m_variations;
  p4compv_t m_eta, m_phi, m_mass;
};

class ModifiedPtMCollection { // map of variation collections
public:
  using LorentzVector = ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float>>;
  using p4compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMCollection() = default;
  ModifiedPtMCollection(p4compv_t eta, p4compv_t phi) : m_eta(eta), m_phi(phi) {}
  void add(const std::string& name, p4compv_t pt, p4compv_t mass)
  {
    m_variations[name] = std::make_pair(pt, mass);
  }

  const p4compv_t& eta() const { return m_eta; }
  const p4compv_t& phi() const { return m_phi; }

  class Variation { // single variation
  public:
    Variation(const ModifiedPtMCollection& parent, const p4compv_t& pt, const p4compv_t& mass) : m_parent(parent), m_pt(pt), m_mass(mass) {}
    const p4compv_t& pt  () const { return m_pt; }
    const p4compv_t& eta () const { return m_parent.eta(); }
    const p4compv_t& phi () const { return m_parent.phi(); }
    const p4compv_t& mass() const { return m_mass; }
    std::size_t size() const { return m_pt.size(); }
  private:
    const ModifiedPtMCollection& m_parent;
    const p4compv_t& m_pt;
    const p4compv_t& m_mass;
  };

  Variation at(const std::string& name) const {
    const auto& varPM = m_variations.at(name);
    return Variation(*this, varPM.first, varPM.second);
  }
  std::size_t size() const { return m_variations.size(); }
private:
  boost::container::flat_map<std::string,std::pair<p4compv_t,p4compv_t>> m_variations;
  p4compv_t m_eta, m_phi;
};

class ModifiedMET {
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedMET() = default;
  // initialize with the nominal value for all variations
  ModifiedMET(std::size_t n, float px_nom, float py_nom)
    : m_px(n, px_nom), m_py(n, py_nom) {}

  std::size_t size() const { return m_px.size(); }
  const compv_t& px() const { return m_px; }
  const compv_t& py() const { return m_py; }
  float px (std::size_t i) const { return m_px[i]; }
  float py (std::size_t i) const { return m_py[i]; }
  float pt (std::size_t i) const { return std::sqrt(m_px[i]*m_px[i]+m_py[i]*m_py[i]); }
  float phi(std::size_t i) const { return std::atan2(m_py[i], m_px[i]); }

  void setXY(std::size_t i, float dpx, float dpy) {
    m_px[i] = dpx;
    m_py[i] = dpy;
  }
  void addR_proj(std::size_t i, float cosphi, float sinphi, float dp) {
    m_px[i] += dp*cosphi;
    m_py[i] += dp*sinphi;
  }
private:
  compv_t m_px;
  compv_t m_py;
};
};
