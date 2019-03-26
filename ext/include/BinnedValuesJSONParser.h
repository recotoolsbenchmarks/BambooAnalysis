#pragma once

class BinnedValues;

namespace BinnedValuesJSONParser {
    BinnedValues parse_file(const std::string& file);
};
