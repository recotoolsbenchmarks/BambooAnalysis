#pragma once

#include <memory>
#include <string>
#include <vector>

namespace torch {
  namespace jit {
    namespace script {
      class Module;
    }
  }
}

namespace bamboo {
class TorchEvaluator {
public:
  using input_t = std::vector<float>;
  using output_t = std::vector<float>;

  //disabled to avoid std::string in the ABI, until pytorch is also built with the C++11 ABI
  //explicit TorchEvaluator( const std::string& scriptName );
  explicit TorchEvaluator( const char* scriptName ); // workaround
  ~TorchEvaluator() noexcept;

  output_t operator() ( const input_t& input ) const;
  output_t operator() ( input_t&& input ) const;
private:
  std::unique_ptr<torch::jit::script::Module> m_module;
};
}
