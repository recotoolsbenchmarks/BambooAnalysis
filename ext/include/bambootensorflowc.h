#pragma once

#include <string>
#include <vector>

extern "C" {
class TF_Graph;
class TF_Session;
class TF_SessionOptions;
class TF_Operation;
typedef struct TF_Output_ {
  TF_Operation* oper;
  int index;  // The index of the output within oper.
} TF_Output_;
class TF_Tensor;
}

namespace bamboo {
class TensorflowCEvaluator {
public:
  using input_t = std::vector<float>;
  using output_t = std::vector<float>;

  explicit TensorflowCEvaluator( const std::string& modelFile,
      const std::string& inputName, const std::string& outputName,
      int64_t nIn, int64_t nOut);
  ~TensorflowCEvaluator();

  output_t operator() ( const input_t& input ) const;
private:
  TF_Graph* m_graph;
  TF_Session* m_session;
  TF_SessionOptions* m_sessionOpts;
  TF_Output_ m_inputs[1];
  TF_Output_ m_outputs[1];
  TF_Tensor* m_inputValues[1];
  int64_t m_nIn, m_nOut;
};
}
