#include "bamboorandom.h"

#include "TRandom3.h"

TRandom3& rdfhelpers::getTRandom3(uint32_t seed) {
  static thread_local TRandom3 rg{};
  rg.SetSeed(seed);
  return rg;
}

std::mt19937& rdfhelpers::getStdMT19937(uint32_t seed) {
  static thread_local std::mt19937 rg{};
  rg.seed(seed);
  return rg;
}
