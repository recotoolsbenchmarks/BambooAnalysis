#include "bamboolwtnn.h"

#include <fstream>
#include "lwtnn/LightweightGraph.hh"
#include "lwtnn/parse_json.hh"

bamboo::LwtnnEvaluator::LwtnnEvaluator(const std::string& fileName, const std::vector<std::pair<std::string,std::string>>& inputNames, const std::vector<std::string>& outputNames)
  : m_inputNames(inputNames), m_outputNames(outputNames)
{
  std::ifstream input(fileName);
  m_lwGraph = std::make_unique<lwt::LightweightGraph>(lwt::parse_json_graph(input));
}

bamboo::LwtnnEvaluator::~LwtnnEvaluator() noexcept = default;

bamboo::LwtnnEvaluator::output_t bamboo::LwtnnEvaluator::evaluate( const input_t& input ) const
{
  if ( input.size() != m_inputNames.size() ) {
    throw std::runtime_error("Input sizes do not match");
  }
  std::map<std::string,std::map<std::string,double>> inputValues;
  for ( std::size_t i{0}; i != m_inputNames.size(); ++i ) {
    const auto& ndNm = m_inputNames[i];
    inputValues[ndNm.first][ndNm.second] = input[i];
  }
  auto outputValues = m_lwGraph->compute(inputValues);
  output_t output;
  for ( const auto& outName : m_outputNames ) {
    output.push_back(outputValues[outName]);
  }
  return output;
}
