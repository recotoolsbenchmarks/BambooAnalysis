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

  output_t evaluate( input_t&& input ) const;
  output_t evaluate( const input_t& input ) const
  {
    input_t input_copy{input};
    return evaluate(std::move(input_copy));
  }

  template<typename RANGE>
  output_t evaluate( RANGE range ) const {
    input_t input;
    std::copy(std::begin(range), std::end(range), std::back_inserter(input));
    return evaluate(std::move(input));
  }
private:
  std::unique_ptr<torch::jit::script::Module> m_module;
};
}
