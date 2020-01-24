#pragma once

#include <random>
// forward declarations
class TRandom3;

/**
 * Thread-local (thread-safe and efficient) random generators
 */
namespace rdfhelpers {
  TRandom3& getTRandom3(uint32_t seed);
  std::mt19937& getStdMT19937(uint32_t seed);
};
