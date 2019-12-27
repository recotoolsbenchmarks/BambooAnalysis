#pragma once

#include <ROOT/RVec.hxx>

namespace rdfhelpers {

class ModifiedPtCollection {
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtCollection() = default;
  ModifiedPtCollection(std::size_t n, const compv_t& pt) : m_pt(n, pt) {}

  std::size_t size() const { return m_pt.size(); }
  const compv_t& pt(std::size_t i) const { return m_pt[i]; }

  void set(std::size_t i, const compv_t& pt) { m_pt[i] = pt; }
  void set(std::size_t i, compv_t&& pt) { m_pt[i] = std::move(pt); }
private:
  std::vector<compv_t> m_pt;
};

class ModifiedPtMCollection { // map of variation collections
public:
  using compv_t = ROOT::VecOps::RVec<float>;

  ModifiedPtMCollection() = default;
  ModifiedPtMCollection(std::size_t n, const compv_t& pt, const compv_t& mass)
    : m_pt(n, pt), m_mass(n, mass) {}

  std::size_t size() const { return m_pt.size(); }

  const compv_t& pt(std::size_t i) const { return m_pt[i]; }
  const compv_t& mass(std::size_t i) const { return m_mass[i]; }

  void set(std::size_t i, const compv_t& pt, const compv_t& mass) {
    m_pt[i] = pt;
    m_mass[i] = mass;
  }
  void set(std::size_t i, compv_t&& pt, compv_t&& mass) {
    m_pt[i] = std::move(pt);
    m_mass[i] = std::move(mass);
  }
private:
  std::vector<compv_t> m_pt;
  std::vector<compv_t> m_mass;
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
