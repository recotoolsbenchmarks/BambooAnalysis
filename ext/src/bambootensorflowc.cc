#include "bambootensorflowc.h"

#include <iostream>
#include <fstream>
#include <sstream>
#include <cstring>
#include <numeric>

namespace {
std::vector<int64_t> getInOutDimensions(TF_Graph* graph, TF_Output node, const std::string& ndName, TF_Status* status) {
  const auto nDims = TF_GraphGetTensorNumDims(graph, node, status);
  if (TF_GetCode(status) != TF_OK) {
    std::stringstream msg;
    msg << "ERROR: could not get the number of dimensions for node " << ndName << " " << TF_Message(status);
    throw std::runtime_error(msg.str());
  }
  std::vector<int64_t> dims(std::size_t(nDims), 0);
  TF_GraphGetTensorShape(graph, node, dims.data(), nDims, status);
  if ((TF_GetCode(status) != TF_OK) || dims.empty() || (dims[0] != -1) ) {
    std::stringstream msg;
    if (TF_GetCode(status) != TF_OK) {
      msg << "ERROR: could not get dimensions for node " << ndName << " " << TF_Message(status);
    } else {
      msg << "ERROR: dimensions for node " << ndName << " are unexpected: (";
      for ( auto dim : dims ) {
        msg << dim << ", ";
      }
      msg << ")";
    }
    throw std::runtime_error(msg.str());
  }
  return dims;
}
void printShape(const std::string& nodeName, const std::vector<int64_t> dims) {
  std::cout << "Tensor shape for node " << nodeName << ": (";
  for ( auto dm : dims ) {
    std::cout << dm << ", ";
  }
  std::cout << ")" << std::endl;
}
}

bamboo::TensorflowCEvaluator::TensorflowCEvaluator(const std::string& modelFile,
    const std::vector<std::string>& inputNames, const std::vector<std::string>& outputNames)
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
  m_inputs.reserve(inputNames.size());
  m_inputValues.reserve(inputNames.size());
  for ( const auto& ndName : inputNames ) {
    auto node = TF_GraphOperationByName(m_graph, ndName.c_str());
    if ( ! node ) {
      throw std::runtime_error("Could not retrieve input node "+ndName);
    }
    m_inputs.push_back({ node, 0 });
    auto dims = getInOutDimensions(m_graph, m_inputs.back(), ndName, status);
    //printShape(ndName, dims);
    dims[0] = 1;
    const auto inSize = std::accumulate(std::begin(dims)+1, std::end(dims), 1,
            [] ( std::size_t prod, int64_t dimSz ) { return prod*dimSz; });
    m_inputValues.push_back(TF_AllocateTensor(TF_FLOAT, dims.data(), dims.size(), sizeof(float)*inSize));
    m_inSize.push_back(inSize);
  }
  m_outputs.reserve(outputNames.size());
  for ( const auto& ndName : outputNames ) {
    auto node = TF_GraphOperationByName(m_graph, ndName.c_str());
    if ( ! node ) {
      throw std::runtime_error("Could not retrieve output node "+ndName);
    }
    m_outputs.push_back({ node, 0 });
    //printShape(ndName, getInOutDimensions(m_graph, m_outputs.back(), ndName, status));
  }
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
  auto itBegin = input.cbegin();
  for ( std::size_t i = 0; i != m_inputs.size(); ++i ) {
    const auto itEnd = itBegin+m_inSize[i];
    if ( ( itEnd == input.cend() ) != ( i+1 == m_inputs.size() ) ) {
      std::stringstream msg;
      msg << "Incorrect number of input size: expected ";
      std::size_t tot = 0;
      for ( auto ndSz : m_inSize ) {
        msg << (tot!=0 ? "+" : "") << std::to_string(ndSz);
        tot += ndSz;
      }
      msg << "=" << std::to_string(tot) << " values, but received " << std::to_string(input.size());
      throw std::runtime_error(msg.str());
    }
    std::copy(itBegin, itEnd, static_cast<float*>(TF_TensorData(m_inputValues[i])));
    itBegin = itEnd;
  }
  auto output_values = std::vector<TF_Tensor*>(m_outputs.size(), nullptr);
  TF_SessionRun(
      m_session,
      nullptr, // run_options
      m_inputs.data(), m_inputValues.data(), m_inputs.size(),
      m_outputs.data(), output_values.data(), m_outputs.size(),
      nullptr, 0, // target opers, ntargets
      nullptr, // run metadata
      status);
  if (TF_GetCode(status) != TF_OK) {
    throw std::runtime_error(std::string{"Problem in evaluating"}+TF_Message(status));
  }
  std::size_t nOut = 0;
  for ( auto out : output_values ) {
    nOut += TF_TensorElementCount(out);
  }
  bamboo::TensorflowCEvaluator::output_t output(nOut, 0.);
  std::size_t start = 0;
  for ( auto out : output_values ) {
    const auto nElm = TF_TensorElementCount(out);
    std::memcpy(output.data()+start, TF_TensorData(out), nElm*sizeof(float));
    TF_DeleteTensor(out);
    start += nElm;
  }
  TF_DeleteStatus(status);
  return output;
}
