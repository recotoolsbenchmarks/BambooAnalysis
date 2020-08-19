#pragma once

#include <string>
#include <vector>

#include "tensorflow/c/c_api.h"

namespace bamboo {
class TensorflowCEvaluator {
public:
  using input_t = std::vector<float>;
  using output_t = std::vector<float>;

  TensorflowCEvaluator(const std::string& modelFile,
      const std::vector<std::string>& inputNames,
      const std::vector<std::string>& outputNames);
  TensorflowCEvaluator(const std::string& modelFile,
      const std::vector<std::string>& inputNames,
      const std::string& outputName)
    : TensorflowCEvaluator(modelFile, inputNames, std::vector<std::string>({ outputName })) {}
  TensorflowCEvaluator(const std::string& modelFile,
      const std::string& inputName,
      const std::vector<std::string>& outputNames)
    : TensorflowCEvaluator(modelFile, std::vector<std::string>({ inputName }), outputNames) {}
  TensorflowCEvaluator(const std::string& modelFile,
      const std::string& inputName,
      const std::string& outputName)
    : TensorflowCEvaluator(modelFile, std::vector<std::string>({ inputName }), std::vector<std::string>({ outputName })) {}
  ~TensorflowCEvaluator();

  output_t evaluate( const input_t& input ) const;

  template<typename RANGE>
  output_t evaluate( RANGE range ) const {
    input_t input;
    std::copy(std::begin(range), std::end(range), std::back_inserter(input));
    return evaluate(input);
  }
private:
  TF_Graph* m_graph;
  TF_Session* m_session;
  TF_SessionOptions* m_sessionOpts;
  std::vector<TF_Output> m_inputs;
  std::vector<TF_Output> m_outputs;
  std::vector<std::size_t> m_inSize;
  std::vector<TF_Tensor*> m_inputValues;
};
}
