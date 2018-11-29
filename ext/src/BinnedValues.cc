#include "BinnedValues.h"

Parameters::Parameters(Parameters&& rhs) {
    m_values = std::move(rhs.m_values);
}

Parameters::Parameters(std::initializer_list<typename value_type::value_type> init) {
    for (auto& i: init) {
        set(i.first, i.second);
    }
}

Parameters& Parameters::setPt(float pt) {
    m_values[BinningVariable::Pt] = pt;
    return *this;
}

Parameters& Parameters::setEta(float eta) {
    m_values[BinningVariable::Eta] = eta;
    m_values[BinningVariable::AbsEta] = fabs(eta);
    return *this;
}

Parameters& Parameters::setBTagDiscri(float d) {
    m_values[BinningVariable::BTagDiscri] = d;
    return *this;
}

Parameters& Parameters::set(const BinningVariable& bin, float value) {
    m_values.emplace(bin, value);

    // Special case for eta
    if (bin == BinningVariable::Eta) {
        m_values.emplace(BinningVariable::AbsEta, std::abs(value));
    }

    return *this;
}

Parameters& Parameters::set(const typename value_type::value_type& value) {
    set(value.first, value.second);
    return *this;
}

std::vector<float> Parameters::toArray(const std::vector<BinningVariable>& binning) const {
    std::vector<float> values;
    for (const auto& bin: binning) {
        const auto& it = m_values.find(bin);
        if (it == m_values.cend()) {
            std::string message{"Parametrisation depends on '" +
                    BinnedValues::variable_to_string_mapping.left.at(bin) +
                    "' but no value for this parameter has been specified. Please call the appropriate 'set' function of the Parameters object"};
            throw std::invalid_argument(message);
        }

        values.push_back(it->second);
    } 

    return values;
}

const BinnedValues::mapping_bimap BinnedValues::variable_to_string_mapping = {
    {BinningVariable::Pt, "Pt"}, {BinningVariable::Eta, "Eta"},
    {BinningVariable::AbsEta, "AbsEta"}, {BinningVariable::BTagDiscri, "BTagDiscri"}
};

void BinnedValues::setVariables(const std::vector<std::string>& v) {
    binning_variables.clear();

    for (const auto& var: v) {
        auto it = variable_to_string_mapping.right.find(var);
        if (it == variable_to_string_mapping.right.end()) {
            throw std::runtime_error("Unknown kind of variable: " + var);
        }

        binning_variables.push_back(it->second);
    }
}
