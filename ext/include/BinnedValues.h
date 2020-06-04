#pragma once

#include "Histogram.h"

#include <memory>
#include <unordered_map>
#include <sstream>

#include <TFormula.h> // TODO move to FormulaEvaluator, and bundle with jet classes

// FIXME lots of global namespace pollution

enum SystVariation {
    Nominal = 0,
    Down = 1,
    Up = 2
};

template <typename T, typename U>
struct bimap {
    typedef std::unordered_map<T, U> left_type;
    typedef std::unordered_map<U, T> right_type;

    left_type left;
    right_type right;
    
    bimap(std::initializer_list<typename left_type::value_type> l) {
        for (auto& v: l) {
            left.insert(v);
            right.insert(typename right_type::value_type(v.second, v.first));
        }
    }

    bimap() {
        // Empty
    }

    bimap(bimap&& rhs) {
        left = std::move(rhs.left);
        right = std::move(rhs.right);
    }
};

/**
 * Enum of possible Variables used for binning values
 */
enum class BinningVariable {
    Pt,
    Eta,
    AbsEta,
    BTagDiscri,
    NumJets,
    NumTrueInteractions
};

// Hash function for BinningVariable enum class
namespace std {
    template<>
    struct hash<BinningVariable> {
        typedef BinningVariable argument_type;
        typedef std::size_t result_type;

        hash<uint8_t> int_hash;
 
        result_type operator()(argument_type const& s) const {
            return int_hash(static_cast<uint8_t>(s));
        }
    };
};

class Parameters {
    public:
        typedef std::unordered_map<BinningVariable, float> value_type;

        Parameters() = default;
        Parameters(Parameters&& rhs);
        Parameters(std::initializer_list<typename value_type::value_type> init);


        Parameters& setPt(float pt);
        Parameters& setEta(float eta);
        Parameters& setBTagDiscri(float d);
        Parameters& set(const BinningVariable& bin, float value);
        Parameters& set(const typename value_type::value_type& value);

        std::vector<float> toArray(const std::vector<BinningVariable>&) const;

    private:
        value_type m_values;
};

class BinnedValues {

    public:
    typedef bimap<BinningVariable, std::string> mapping_bimap;
    static const mapping_bimap variable_to_string_mapping;

    /**
     * Type of possible errors: Suppose we have E +- ΔE
     *   - ABSOLUTE = ΔE
     *   - RELATIVE = ΔE / E
     *   - VARIATED = E + ΔE or E - ΔE
     */
    enum class ErrorType {
        ABSOLUTE,
        RELATIVE,
        VARIATED
    };

    BinnedValues(BinnedValues&& rhs) = default;

    BinnedValues() = default;

    const Histogram<float>& binned() const { return *m_binned; }
    Histogram<float>& binned() { return *m_binned; }
    void setBinned(std::unique_ptr<Histogram<float>>&& binned) { m_binned = std::move(binned); }

    const Histogram<std::shared_ptr<TFormula>, float>& formula() const { return *m_formula; }
    Histogram<std::shared_ptr<TFormula>, float>& formula() { return *m_formula; }
    void setFormula(std::unique_ptr<Histogram<std::shared_ptr<TFormula>, float>>&& formula) { m_formula = std::move(formula); }

    void setVariables(const std::vector<std::string>&);
    void setFormulaVariableIndex(std::size_t idx)
    {
      use_formula = true;
      formula_variable_index = idx;
    }
    void setRange(float xMin, float xMax)
    {
      minimum = xMin;
      maximum = xMax;
    }
    void setErrorType(ErrorType eType) { error_type = eType; }

    private:
    template <typename _Value>
        std::vector<_Value> get(Histogram<_Value, float>& h, const std::vector<float>& bins, bool& outOfRange) const {
            std::size_t bin = h.findClosestBin(bins, &outOfRange);
            if (bin == 0) {
                std::stringstream msg;
                msg << "Failed to found the right bin for a scale-factor. This should not happend. Bins: [";
                for (float b: bins) {
                    msg << b << ", ";
                }

                msg.seekp(msg.tellp() - 2l);
                msg << "]";

                throw std::runtime_error(msg.str());
            }

            return {h.getBinContent(bin), h.getBinErrorLow(bin), h.getBinErrorHigh(bin)};
        }

    // List of variables used in the binning. First entry is the X variable, second one Y, etc.
    std::vector<BinningVariable> binning_variables;

    bool use_formula = false;

    // Binned data
    std::unique_ptr<Histogram<float>> m_binned;

    // Formula data
    std::unique_ptr<Histogram<std::shared_ptr<TFormula>, float>> m_formula;

    ErrorType error_type;
    size_t formula_variable_index = -1; // Only used in formula mode

    float maximum;
    float minimum;

    /**
     * Convert relative errors to absolute errors
     **/
    std::vector<float> relative_errors_to_absolute(const std::vector<float>& array) const {
        std::vector<float> result(3);
        result[Nominal] = array[Nominal];
        result[Up] = array[Nominal] * array[Up];
        result[Down] = array[Nominal] * array[Down];

        return result;
    };

    /**
     * Convert variated errors to absolute errors
     **/
    std::vector<float> variated_errors_to_absolute(const std::vector<float>& array) const {
        std::vector<float> result(3);
        result[Nominal] = array[Nominal];
        result[Up] = std::abs(array[Up] - array[Nominal]);
        result[Down] = std::abs(array[Nominal] - array[Down]);

        return result;
    };

    std::vector<float> convert_errors(const std::vector<float>& array) const {
        switch (error_type) {
            case ErrorType::ABSOLUTE:
                return array;

            case ErrorType::RELATIVE:
                return relative_errors_to_absolute(array);

            case ErrorType::VARIATED:
                return variated_errors_to_absolute(array);
        }

        throw std::runtime_error("Invalid error type");
    }

    /**
     * Check that the up and down variation are
     * still between the allowed range
     **/
    void clamp(std::vector<float>& array) const {
        if ((array[Nominal] + array[Up]) > maximum) {
            array[Up] = maximum - array[Nominal];
        }

        if ((array[Nominal] - array[Down]) < minimum) {
            array[Down] = -(minimum - array[Nominal]);
        }
    }

    public:
    // returns { E, ΔE+, ΔE- }
    virtual std::vector<float> get(const Parameters& parameters) const {
        static auto double_errors = [](std::vector<float>& values) {
            values[Up] *= 2;
            values[Down] *= 2;
        };

        std::vector<float> variables = parameters.toArray(binning_variables);

        bool outOfRange = false;

        if (!use_formula) {
            if (! m_binned)
                return {0., 0., 0.};

            std::vector<float> values = convert_errors(get<float>(*m_binned, variables, outOfRange));

            if (outOfRange)
                double_errors(values);

            clamp(values);

            return values;
        } else {
            if (! m_formula)
                return {0., 0., 0.};

            std::vector<std::shared_ptr<TFormula>> formulas = get<std::shared_ptr<TFormula>>(*m_formula, variables, outOfRange);
            std::vector<float> values;

            // Ensure variables are not outside the validity range
            variables = m_formula->clamp(variables);

            for (auto& formula: formulas) {
                values.push_back(formula->Eval(variables[formula_variable_index]));
            }

            values = convert_errors(values);

            if (outOfRange)
                double_errors(values);

            clamp(values);

            return values;
        }
    }

};
