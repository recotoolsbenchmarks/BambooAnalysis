#include <boost/property_tree/json_parser.hpp>

#include "BinnedValuesJSONParser.h"

void BinnedValuesJSONParser::parse_file(const std::string& file) {

    boost::property_tree::ptree ptree;
    boost::property_tree::read_json(file, ptree);

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

    m_values.setVariables(variables);

    bool formula = ptree.get<bool>("formula", false);
    
    m_values.use_formula = formula;
    if (formula) {
        std::string variable = ptree.get<std::string>("variable");

        if (variable == "x")
            m_values.formula_variable_index = 0;
        else if (variable == "y")
            m_values.formula_variable_index = 1;
        else if (variable == "z")
            m_values.formula_variable_index = 2;
        else {
            std::string message{"Unsupported variable: " + variable};
            throw std::logic_error(message);
        }
    }

    m_values.maximum = ptree.get("maximum", std::numeric_limits<float>::max());
    m_values.minimum = ptree.get("minimum", 0);

    std::string error_type = ptree.get<std::string>("error_type");
    std::transform(error_type.begin(), error_type.end(), error_type.begin(), ::tolower);

    if (error_type == "absolute")
        m_values.error_type = BinnedValues::ErrorType::ABSOLUTE;
    else if (error_type == "relative")
        m_values.error_type = BinnedValues::ErrorType::RELATIVE;
    else if (error_type == "variated")
        m_values.error_type = BinnedValues::ErrorType::VARIATED;
    else
        throw std::runtime_error("Invalid error_type. Only 'absolute', 'relative' and 'variated' are supported");

    switch (dimension) {
        case 1:
            if (!formula)
                m_values.binned.reset(new OneDimensionHistogram<float>(binning_x));
            else
                m_values.formula.reset(new OneDimensionHistogram<std::shared_ptr<TFormula>, float>(binning_x));
            break;

        case 2:
            if (!formula)
                m_values.binned.reset(new TwoDimensionsHistogram<float>(binning_x, binning_y));
            else
                m_values.formula.reset(new TwoDimensionsHistogram<std::shared_ptr<TFormula>, float>(binning_x, binning_y));
            break;

        case 3:
            if (!formula)
                m_values.binned.reset(new ThreeDimensionsHistogram<float>(binning_x, binning_y, binning_z));
            else
                m_values.formula.reset(new ThreeDimensionsHistogram<std::shared_ptr<TFormula>, float>(binning_x, binning_y, binning_z));
            break;
    }

    if (formula) {
        parse_data<std::string>(ptree, dimension);
    } else {
        parse_data<float>(ptree, dimension);
    }
}
