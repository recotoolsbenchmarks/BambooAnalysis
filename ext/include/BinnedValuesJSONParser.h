#pragma once

#include "Histogram.h"
#include "BinnedValues.h"

#include <boost/property_tree/ptree.hpp>

#include <memory>

#include <TFormula.h>

class BinnedValuesJSONParser {

    public:
        BinnedValuesJSONParser(const std::string& file) {
            parse_file(file);
        }

        virtual BinnedValues&& get_values() final {
            return std::move(m_values);
        }

    private:
        void parse_file(const std::string& file);

        std::vector<float> get_array(boost::property_tree::ptree& ptree) {
            std::vector<float> vector;
            for (auto& value: ptree) {
                vector.push_back(std::stof(value.second.data()));
            }

            return vector;
        }

        std::vector<std::string> get_string_array(boost::property_tree::ptree& ptree) {
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
            fillHistogram(*val.binned.get(), bins, value, error_low, error_high);
        }

        void fillHistogram(BinnedValues& val, const std::vector<float>& bins, const std::string& value, const std::string& error_low, const std::string& error_high) {
            std::shared_ptr<TFormula> value_formula(new TFormula("", value.c_str()));
            std::shared_ptr<TFormula> error_low_formula(new TFormula("", error_low.c_str()));
            std::shared_ptr<TFormula> error_high_formula(new TFormula("", error_high.c_str()));
            fillHistogram(*val.formula.get(), bins, value_formula, error_low_formula, error_high_formula);
        }

        template <typename _Content>
        void parse_data(boost::property_tree::ptree& ptree, std::size_t dimension) {
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

                                fillHistogram(m_values, {mean_x, mean_y, mean_z}, value, error_low, error_high);
                            }

                        } else {

                            _Content value = data_y.second.get<_Content>("value");
                            _Content error_low = data_y.second.get<_Content>("error_low");
                            _Content error_high = data_y.second.get<_Content>("error_high");

                            fillHistogram(m_values, {mean_x, mean_y}, value, error_low, error_high);
                        }

                    }

                } else {

                    _Content value = data_x.second.get<_Content>("value");
                    _Content error_low = data_x.second.get<_Content>("error_low");
                    _Content error_high = data_x.second.get<_Content>("error_high");

                    fillHistogram(m_values, {mean_x}, value, error_low, error_high);
                }
            }
        }

        BinnedValues m_values;
};
