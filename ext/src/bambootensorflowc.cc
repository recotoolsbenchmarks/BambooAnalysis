#include "bambootensorflowc.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include "tensorflow/c/c_api.h"

bamboo::TensorflowCEvaluator::TensorflowCEvaluator(const std::string& modelFile,
    const std::string& inputName, const std::string& outputName,
    int64_t nIn, int64_t nOut )
  : m_nIn(nIn), m_nOut(nOut)
{
  std::ifstream input{modelFile, std::ios::binary};
  std::vector<unsigned char> inBuffer(std::istreambuf_iterator<char>(input), {});
  TF_Buffer* graphDef = TF_NewBuffer();
  graphDef->data = inBuffer.data();
  graphDef->length = inBuffer.size();
  TF_Status* status = TF_NewStatus();
  TF_ImportGraphDefOptions* importOpts = TF_NewImportGraphDefOptions();
  m_graph = TF_NewGraph();
  TF_GraphImportGraphDef(m_graph, graphDef, importOpts, status);
  TF_DeleteImportGraphDefOptions(importOpts);
  if (TF_GetCode(status) != TF_OK) {
    throw std::runtime_error(std::string{"ERROR: could not import graph "}+TF_Message(status));
  }
  m_sessionOpts = TF_NewSessionOptions();
  m_session = TF_NewSession(m_graph, m_sessionOpts, status);
  if (TF_GetCode(status) != TF_OK) {
    throw std::runtime_error(std::string{"ERROR: could not create session "}+TF_Message(status));
  }
  auto inNode = TF_GraphOperationByName(m_graph, inputName.c_str());
  if ( ! inNode ) {
    throw std::runtime_error("Could not retrieve input node "+inputName);
  }
  auto outNode = TF_GraphOperationByName(m_graph, outputName.c_str());
  if ( ! outNode ) {
    throw std::runtime_error("Could not retrieve output node "+outputName);
  }
  m_inputs[0] = { inNode, 0 };
  m_outputs[0] = { outNode, 0 };
  const int64_t inDims[2] = { 1, m_nIn };
  m_inputValues[0] = TF_AllocateTensor(TF_FLOAT, inDims, 2, sizeof(float)*1*m_nIn);
  TF_DeleteStatus(status);
}

bamboo::TensorflowCEvaluator::~TensorflowCEvaluator()
{
  TF_Status* status = TF_NewStatus();
  TF_CloseSession(m_session, status);
  TF_DeleteSession(m_session, status);
  TF_DeleteSessionOptions(m_sessionOpts);
  TF_DeleteGraph(m_graph);
  TF_DeleteStatus(status);
}

bamboo::TensorflowCEvaluator::output_t bamboo::TensorflowCEvaluator::evaluate( const bamboo::TensorflowCEvaluator::input_t& input ) const
{
  TF_Status* status = TF_NewStatus();
  if ( m_nIn != input.size() ) {
    throw std::runtime_error("Incorrect number of input values: expected "+std::to_string(m_nIn)+" but received "+std::to_string(input.size()));
  }
  std::copy(std::begin(input), std::end(input), static_cast<float*>(TF_TensorData(m_inputValues[0])));
  TF_Tensor* output_values[1] = { nullptr };
  TF_SessionRun(
      m_session,
      nullptr, // run_options
      reinterpret_cast<const TF_Output*>(m_inputs), m_inputValues, 1,
      reinterpret_cast<const TF_Output*>(m_outputs), output_values, 1,
      nullptr, 0, // target opers, ntargets
      nullptr, // run metadata
      status);
  if (TF_GetCode(status) != TF_OK) {
    throw std::runtime_error(std::string{"Problem in evaluating"}+TF_Message(status));
  }
  if ( ( 2 != TF_NumDims(output_values[0]) ) || ( 1 != TF_Dim(output_values[0], 0) ) || ( m_nOut != TF_Dim(output_values[0], 1) ) ) {
    std::stringstream outShape{};
    outShape << TF_Dim(output_values[0], 0);
    for ( std::size_t i{1}; i != TF_NumDims(output_values[0]); ++i ) {
      outShape << "x" << TF_Dim(output_values[0], i);
    }
    throw std::runtime_error("Incorrect output shape: expected 1x"+std::to_string(m_nOut)+" but received "+outShape.str());
  }
  bamboo::TensorflowCEvaluator::output_t output(std::size_t(m_nOut), 0.);
  std::memcpy(output.data(), TF_TensorData(output_values[0]), m_nOut*sizeof(float));
  TF_DeleteTensor(output_values[0]);
  TF_DeleteStatus(status);
  return output;
}
