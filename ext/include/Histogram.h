#pragma once

#include <memory>
#include <vector>
#include <algorithm>

template<typename T, typename _Bin = T>
class Histogram {
    public:
        virtual ~Histogram();

        virtual std::size_t findBin(const std::vector<_Bin>& values) = 0;
        virtual std::size_t findClosestBin(const std::vector<_Bin>& values, bool* outOfRange = nullptr) = 0;
        virtual bool inRange(const std::vector<_Bin>& values) = 0;
        virtual std::vector<_Bin> clamp(const std::vector<_Bin>& values) = 0;

        T getBinContent(std::size_t bin) {
            return m_values[bin - 1];
        }
        T getBinErrorLow(std::size_t bin) {
            return m_errors_low[bin - 1];
        }
        T getBinErrorHigh(std::size_t bin) {
            return m_errors_high[bin - 1];
        }

        void setBinContent(std::size_t bin, T value) {
            m_values[bin - 1] = value;
        }
        void setBinErrorLow(std::size_t bin, T value) {
            m_errors_low[bin - 1] = value;
        }
        void setBinErrorHigh(std::size_t bin, T value) {
            m_errors_high[bin - 1] = value;
        }

        void setContent(const std::vector<_Bin>& values, T content) {
            size_t bin = findBin(values);
            if (bin == 0)
                return;

            setBinContent(bin, content);
        }

        size_t size() const {
            return m_size;
        }

    protected:
        Histogram(std::size_t size) {
            m_size = size;

            m_values.reset(new T[m_size]);
            m_errors_low.reset(new T[m_size]);
            m_errors_high.reset(new T[m_size]);
        }

        static size_t findBin(const std::vector<_Bin>& array, _Bin value) {
            return ( ( value < array[0] ) || ( value >= array[array.size()-1] ) ) ? 0 :
                std::distance(array.begin(), std::upper_bound(array.begin(), array.end(), value));
        }

        static size_t findClosestBin(const std::vector<_Bin>& array, _Bin value, bool* outOfRange = nullptr) {
            if (outOfRange)
                *outOfRange = false;

            if (value < array.front()) {
                if (outOfRange)
                    *outOfRange = true;
                return 1;
            } else if (value >= array.back()) {
                if (outOfRange)
                    *outOfRange = true;
                return array.size() - 1;
            } else {
                return findBin(array, value);
            }
        }

        static bool inRange(const std::vector<_Bin>& array, _Bin value) {
            _Bin min = array.front();
            _Bin max = array.back();

            return ((value >= min) && (value < max));
        }

        static _Bin clamp(const std::vector<_Bin>& array, _Bin value) {
            _Bin min = array.front();
            _Bin max = array.back();

            if (value < min)
                return min;

            if (value > max)
                return max;

            return value;
        }

        std::size_t m_size;
        std::unique_ptr<T[]> m_values;
        std::unique_ptr<T[]> m_errors_low;
        std::unique_ptr<T[]> m_errors_high;

    private:
        Histogram() = delete;

};

template<typename T, typename _Bin = T>
class OneDimensionHistogram: public Histogram<T, _Bin> {
    public:
        OneDimensionHistogram(const std::vector<_Bin>& bins):
            Histogram<T, _Bin>(bins.size() - 1) {
                m_bins = bins;
        }
        virtual ~OneDimensionHistogram();

        virtual std::size_t findBin(const std::vector<_Bin>& values) override {
            if (values.size() != 1)
                return 0;

            _Bin value = values[0];

            return Histogram<T, _Bin>::findBin(m_bins, value);
        }

        virtual std::size_t findClosestBin(const std::vector<_Bin>& values, bool* outOfRange = nullptr) override {
            if (values.size() != 1)
                return 0;

            _Bin value = values[0];

            return Histogram<T, _Bin>::findClosestBin(m_bins, value, outOfRange);
        }

        virtual bool inRange(const std::vector<_Bin>& values) override {
            if (values.size() != 1) {
                return false;
            }

            _Bin value = values.front();
            return Histogram<T, _Bin>::inRange(m_bins, value);
        }

        virtual std::vector<_Bin> clamp(const std::vector<_Bin>& values) override {
            if (values.size() != 1) {
                return values;
            }

            _Bin value = values.front();
            return {Histogram<T, _Bin>::clamp(m_bins, value)};
        }

    private:
        std::vector<_Bin> m_bins;
};

template<typename T, typename _Bin = T>
class TwoDimensionsHistogram: public Histogram<T, _Bin> {
    public:
        TwoDimensionsHistogram(const std::vector<_Bin>& bins_x, const std::vector<_Bin>& bins_y):
            Histogram<T, _Bin>((bins_x.size() - 1) * (bins_y.size() - 1)) {
                m_bins_x = bins_x;
                m_bins_y = bins_y;
        }
        virtual ~TwoDimensionsHistogram();

        virtual std::size_t findBin(const std::vector<_Bin>& values) override {
            if (values.size() != 2)
                return 0;

            _Bin value_x = values[0];
            _Bin value_y = values[1];

            size_t bin_x = Histogram<T, _Bin>::findBin(m_bins_x, value_x);
            if (bin_x == 0)
                return 0;

            size_t bin_y = Histogram<T, _Bin>::findBin(m_bins_y, value_y);
            if (bin_y == 0)
                return 0;

            return bin_x + (m_bins_x.size() - 1) * (bin_y - 1);
        }

        virtual std::size_t findClosestBin(const std::vector<_Bin>& values, bool* outOfRange = nullptr) override {
            if (values.size() != 2)
                return 0;

            _Bin value_x = values[0];
            _Bin value_y = values[1];
            bool local_outOfRange = false;

            if (outOfRange)
                *outOfRange = false;

            size_t bin_x = Histogram<T, _Bin>::findClosestBin(m_bins_x, value_x, &local_outOfRange);

            if (outOfRange)
                *outOfRange |= local_outOfRange;

            if (bin_x == 0)
                return 0;

            size_t bin_y = Histogram<T, _Bin>::findClosestBin(m_bins_y, value_y, &local_outOfRange);


            if (outOfRange)
                *outOfRange |= local_outOfRange;

            if (bin_y == 0)
                return 0;

            return bin_x + (m_bins_x.size() - 1) * (bin_y - 1);
        }

        virtual bool inRange(const std::vector<_Bin>& values) override {
            if (values.size() != 2) {
                return false;
            }

            _Bin value_x = values.front();
            _Bin value_y = values[1];

            return Histogram<T, _Bin>::inRange(m_bins_x, value_x) && Histogram<T, _Bin>::inRange(m_bins_y, value_y);
        }

        virtual std::vector<_Bin> clamp(const std::vector<_Bin>& values) override {
            if (values.size() != 2) {
                return values;
            }

            _Bin value_x = values.front();
            _Bin value_y = values[1];

            return {Histogram<T, _Bin>::clamp(m_bins_x, value_x), Histogram<T, _Bin>::clamp(m_bins_y, value_y)};
        }

    private:
        std::vector<_Bin> m_bins_x;
        std::vector<_Bin> m_bins_y;
};

template<typename T, typename _Bin = T>
class ThreeDimensionsHistogram: public Histogram<T, _Bin> {
    public:
        ThreeDimensionsHistogram(const std::vector<_Bin>& bins_x, const std::vector<_Bin>& bins_y, const std::vector<_Bin>& bins_z):
            Histogram<T, _Bin>((bins_x.size() - 1) * (bins_y.size() - 1) * (bins_z.size() - 1)) {
                m_bins_x = bins_x;
                m_bins_y = bins_y;
                m_bins_z = bins_z;
        }
        virtual ~ThreeDimensionsHistogram();

        virtual std::size_t findBin(const std::vector<_Bin>& values) override {
            if (values.size() != 3)
                return 0;

            _Bin value_x = values[0];
            _Bin value_y = values[1];
            _Bin value_z = values[2];

            size_t bin_x = Histogram<T, _Bin>::findBin(m_bins_x, value_x);
            if (bin_x == 0)
                return 0;

            size_t bin_y = Histogram<T, _Bin>::findBin(m_bins_y, value_y);
            if (bin_y == 0)
                return 0;

            size_t bin_z = Histogram<T, _Bin>::findBin(m_bins_z, value_z);
            if (bin_z == 0)
                return 0;

            return bin_x + (m_bins_x.size() - 1) * ((bin_y - 1) + (m_bins_y.size() - 1) * (bin_z - 1));
        }

        virtual std::size_t findClosestBin(const std::vector<_Bin>& values, bool* outOfRange = nullptr) override {
            if (values.size() != 3)
                return 0;

            _Bin value_x = values[0];
            _Bin value_y = values[1];
            _Bin value_z = values[2];
            bool local_outOfRange = false;

            size_t bin_x = Histogram<T, _Bin>::findClosestBin(m_bins_x, value_x, &local_outOfRange);

            if (outOfRange)
                *outOfRange |= local_outOfRange;

            if (bin_x == 0)
                return 0;

            size_t bin_y = Histogram<T, _Bin>::findClosestBin(m_bins_y, value_y, &local_outOfRange);

            if (outOfRange)
                *outOfRange |= local_outOfRange;

            if (bin_y == 0)
                return 0;

            size_t bin_z = Histogram<T, _Bin>::findClosestBin(m_bins_z, value_z, &local_outOfRange);

            if (outOfRange)
                *outOfRange |= local_outOfRange;

            if (bin_z == 0)
                return 0;

            return bin_x + (m_bins_x.size() - 1) * ((bin_y - 1) + (m_bins_y.size() - 1) * (bin_z - 1));
        }

        virtual bool inRange(const std::vector<_Bin>& values) override {
            if (values.size() != 3) {
                return false;
            }

            _Bin value_x = values[0];
            _Bin value_y = values[1];
            _Bin value_z = values[2];

            return Histogram<T, _Bin>::inRange(m_bins_x, value_x) && Histogram<T, _Bin>::inRange(m_bins_y, value_y) && Histogram<T, _Bin>::inRange(m_bins_z, value_z);
        }

        virtual std::vector<_Bin> clamp(const std::vector<_Bin>& values) override {
            if (values.size() != 3) {
                return values;
            }

            _Bin value_x = values[0];
            _Bin value_y = values[1];
            _Bin value_z = values[2];

            return {Histogram<T, _Bin>::clamp(m_bins_x, value_x), Histogram<T, _Bin>::clamp(m_bins_y, value_y), Histogram<T, _Bin>::clamp(m_bins_z, value_z)};
        }

    private:
        std::vector<_Bin> m_bins_x;
        std::vector<_Bin> m_bins_y;
        std::vector<_Bin> m_bins_z;
};

template <typename _Value, typename _Bin = _Value>
std::size_t findBin(const Histogram<_Value, _Bin>& object, const std::vector<_Bin>& values) {
    return object.findBin(values);
}

template <typename _Value, typename _Bin = _Value>
_Value getBinContent(const Histogram<_Value, _Bin>& object, std::size_t bin) {
    return object.getBinContent(bin);
}

template <typename _Value, typename _Bin = _Value>
_Value getBinErrorLow(const Histogram<_Value, _Bin>& object, std::size_t bin) {
    return object.getBinErrorLow(bin);
}

template <typename _Value, typename _Bin = _Value>
_Value getBinErrorHigh(const Histogram<_Value, _Bin>& object, std::size_t bin) {
    return object.getBinErrorHigh(bin);
}

template<typename T, typename _Bin> Histogram<T,_Bin>::~Histogram() {}
template<typename T, typename _Bin> OneDimensionHistogram<T,_Bin>::~OneDimensionHistogram() {}
template<typename T, typename _Bin> TwoDimensionsHistogram<T,_Bin>::~TwoDimensionsHistogram() {}
template<typename T, typename _Bin> ThreeDimensionsHistogram<T,_Bin>::~ThreeDimensionsHistogram() {}
