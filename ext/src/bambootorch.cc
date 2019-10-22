#include "bambootorch.h"

#include "torch/script.h"

// disabled to avoid std::string in the ABI, until pytorch is also built with the C++11 ABI - workaround below
//bamboo::TorchEvaluator::TorchEvaluator( const std::string& scriptName )
//  : m_module{std::make_unique<torch::jit::script::Module>(torch::jit::load(scriptName))}
//{}

bamboo::TorchEvaluator::TorchEvaluator( const char* scriptName )
  : m_module{std::make_unique<torch::jit::script::Module>(torch::jit::load(std::string{scriptName}))}
{}

bamboo::TorchEvaluator::~TorchEvaluator() noexcept = default;

bamboo::TorchEvaluator::output_t bamboo::TorchEvaluator::operator() ( bamboo::TorchEvaluator::input_t&& input ) const
{
  std::vector<torch::jit::IValue> inputs{};
  inputs.push_back(torch::from_blob(input.data(), {1,int(input.size())}));

  auto output = m_module->forward(inputs).toTensor();

  auto output_a = output.accessor<float,2>();

  bamboo::TorchEvaluator::output_t out;
  out.reserve(output_a.size(1));
  for ( std::size_t i{0}; i != output_a.size(1); ++i ) {
    out.push_back(output_a[0][i]);
  }
  return out;
}


bamboo::TorchEvaluator::output_t bamboo::TorchEvaluator::operator() ( const bamboo::TorchEvaluator::input_t& input ) const
{
  auto input_copy = input;
  return operator() (std::move(input_copy));
}
