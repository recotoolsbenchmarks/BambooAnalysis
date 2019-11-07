#pragma once

#include <memory>
#include <string>
#include <vector>

namespace lwt {
class LightweightGraph;
}

namespace bamboo {
class LwtnnEvaluator {
public:
  using input_t = std::vector<float>;
  using output_t = std::vector<float>;

  LwtnnEvaluator(const std::string& fileName, const std::vector<std::pair<std::string,std::string>>& inputNames, const std::vector<std::string>& outputNames);
  ~LwtnnEvaluator() noexcept;

  output_t evaluate( const input_t& input ) const;
private:
  std::unique_ptr<lwt::LightweightGraph> m_lwGraph;
  std::vector<std::pair<std::string,std::string>> m_inputNames;
  std::vector<std::string> m_outputNames;
};
}
