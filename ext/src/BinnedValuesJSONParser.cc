#include <memory>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>

#include <TFormula.h>

#include "Histogram.h"
#include "BinnedValues.h"
#include "BinnedValuesJSONParser.h"

namespace {
    std::vector<float> get_array(const boost::property_tree::ptree& ptree) {
        std::vector<float> vector;
        for (auto& value: ptree) {
            vector.push_back(std::stof(value.second.data()));
        }

        return vector;
    }

    std::vector<std::string> get_string_array(const boost::property_tree::ptree& ptree) {
        std::vector<std::string> vector;
        for (auto& value: ptree) {
            vector.push_back(value.second.data());
        }

        return vector;
    }

    template <class T, typename _Value>
    void fillHistogram(T& h, const std::vector<float>& bins, const _Value& value, const _Value& error_low, const _Value& error_high) {
        std::size_t bin = h.findBin(bins);

        h.setBinContent(bin, value);
        h.setBinErrorLow(bin, error_low);
        h.setBinErrorHigh(bin, error_high);
    }

    void fillHistogram(BinnedValues& val, const std::vector<float>& bins, const float& value, const float& error_low, const float& error_high) {
        fillHistogram(val.binned(), bins, value, error_low, error_high);
    }

    void fillHistogram(BinnedValues& val, const std::vector<float>& bins, const std::string& value, const std::string& error_low, const std::string& error_high) {
        std::shared_ptr<TFormula> value_formula(new TFormula("", value.c_str()));
        std::shared_ptr<TFormula> error_low_formula(new TFormula("", error_low.c_str()));
        std::shared_ptr<TFormula> error_high_formula(new TFormula("", error_high.c_str()));
        fillHistogram(val.formula(), bins, value_formula, error_low_formula, error_high_formula);
    }

    template <typename _Content>
    void parse_data(BinnedValues& values, const boost::property_tree::ptree& ptree, std::size_t dimension) {
        for (auto& data_x: ptree.get_child("data")) {
            std::vector<float> binning_x = get_array(data_x.second.get_child("bin"));
            float mean_x = (binning_x[0] + binning_x[1]) / 2.;

            if (dimension > 1) {

                for (auto& data_y: data_x.second.get_child("values")) {
                    std::vector<float> binning_y = get_array(data_y.second.get_child("bin"));
                    float mean_y = (binning_y[0] + binning_y[1]) / 2.;

                    if (dimension > 2) {

                        for (auto& data_z: data_y.second.get_child("values")) {
                            std::vector<float> binning_z = get_array(data_z.second.get_child("bin"));
                            float mean_z = (binning_z[0] + binning_z[1]) / 2.;
                            _Content value = data_z.second.get<_Content>("value");
                            _Content error_low = data_z.second.get<_Content>("error_low");
                            _Content error_high = data_z.second.get<_Content>("error_high");

                            fillHistogram(values, {mean_x, mean_y, mean_z}, value, error_low, error_high);
                        }

                    } else {

                        _Content value = data_y.second.get<_Content>("value");
                        _Content error_low = data_y.second.get<_Content>("error_low");
                        _Content error_high = data_y.second.get<_Content>("error_high");

                        fillHistogram(values, {mean_x, mean_y}, value, error_low, error_high);
                    }

                }

            } else {

                _Content value = data_x.second.get<_Content>("value");
                _Content error_low = data_x.second.get<_Content>("error_low");
                _Content error_high = data_x.second.get<_Content>("error_high");

                fillHistogram(values, {mean_x}, value, error_low, error_high);
            }
        }
    }
}

BinnedValues BinnedValuesJSONParser::parse_file(const std::string& file)
{
    boost::property_tree::ptree ptree;
    boost::property_tree::read_json(file, ptree);

    BinnedValues values{};

    size_t dimension = ptree.get<size_t>("dimension", 1);

    std::vector<std::string> variables = get_string_array(ptree.get_child("variables"));
    if (variables.size() != dimension) {
        std::string message{"Invalid number of variables. Expected " + std::to_string(dimension) + ", got " + std::to_string(variables.size())};
        throw std::logic_error(message);
    }

    std::vector<float> binning_x = get_array(ptree.get_child("binning.x"));

    std::vector<float> binning_y;
    std::vector<float> binning_z;
    if (dimension > 1)
        binning_y = get_array(ptree.get_child("binning.y"));
    if (dimension > 2)
        binning_z = get_array(ptree.get_child("binning.z"));

    values.setVariables(variables);

    const bool formula = ptree.get<bool>("formula", false);
    if ( formula ) {
        std::string variable = ptree.get<std::string>("variable");

        if (variable == "x")
            values.setFormulaVariableIndex(0);
        else if (variable == "y")
            values.setFormulaVariableIndex(1);
        else if (variable == "z")
            values.setFormulaVariableIndex(2);
        else {
            std::string message{"Unsupported variable: " + variable};
            throw std::logic_error(message);
        }
    }

    values.setRange(
        ptree.get("maximum", std::numeric_limits<float>::max()),
        ptree.get("minimum", 0)
        );

    std::string error_type = ptree.get<std::string>("error_type");
    std::transform(error_type.begin(), error_type.end(), error_type.begin(), ::tolower);

    if (error_type == "absolute")
        values.setErrorType(BinnedValues::ErrorType::ABSOLUTE);
    else if (error_type == "relative")
        values.setErrorType(BinnedValues::ErrorType::RELATIVE);
    else if (error_type == "variated")
        values.setErrorType(BinnedValues::ErrorType::VARIATED);
    else
        throw std::runtime_error("Invalid error_type. Only 'absolute', 'relative' and 'variated' are supported");

    switch (dimension) {
        case 1:
            if (!formula)
                values.setBinned(std::make_unique<OneDimensionHistogram<float>>(binning_x));
            else
                values.setFormula(std::make_unique<OneDimensionHistogram<std::shared_ptr<TFormula>, float>>(binning_x));
            break;

        case 2:
            if (!formula)
                values.setBinned(std::make_unique<TwoDimensionsHistogram<float>>(binning_x, binning_y));
            else
                values.setFormula(std::make_unique<TwoDimensionsHistogram<std::shared_ptr<TFormula>, float>>(binning_x, binning_y));
            break;

        case 3:
            if (!formula)
                values.setBinned(std::make_unique<ThreeDimensionsHistogram<float>>(binning_x, binning_y, binning_z));
            else
                values.setFormula(std::make_unique<ThreeDimensionsHistogram<std::shared_ptr<TFormula>, float>>(binning_x, binning_y, binning_z));
            break;
    }

    if (formula) {
        parse_data<std::string>(values, ptree, dimension);
    } else {
        parse_data<float>(values, ptree, dimension);
    }

    return values;
}
