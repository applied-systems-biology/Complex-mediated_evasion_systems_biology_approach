// Copyright by Yann Bachelot
//
// Research Group Applied Systems Biology - Head: Prof. Dr. Marc Thilo Figge
// https://www.leibniz-hki.de/en/applied-systems-biology.html
// HKI-Center for Systems Biology of Infection
// Leibniz Institute for Natural Product Research and Infection Biology - Hans Knöll Institute (HKI)
// Adolf-Reichwein-Straße 23, 07745 Jena, Germany
//
// The project code is licensed under BSD2-Clause.
// See the LICENSE file provided with the code for the full license.

#ifndef CORE_ANALYSER_HISTOGRAM_MEASURMENT_H_
#define CORE_ANALYSER_HISTOGRAM_MEASURMENT_H_

#include <map>
#include <vector>
#include <variant>

#include "core/utils/misc_util.h"

class HistogramMeasurement {
public:
  // Class for  histogram measurements. Histogram measurements are agglomerating information over time.
    using cache_type = std::variant<std::monostate, int, double, std::string>;
    template<typename... Args>
    explicit HistogramMeasurement(std::string id, Args &&... args): owner_id_(std::move(id)),
                                                                    number_of_elements_(sizeof...(Args)),
                                                                    keys_({std::forward<Args>(args)...}) {}

    template<typename... Args>
    typename std::enable_if<(std::is_convertible<Args, std::pair<const char *, cache_type>>::value && ...), void>::type
    addValues(Args &&...args) {
        ((data_[args.first].emplace_back(args.second)), ...);
    }

    template<typename... Args>
    typename std::enable_if<std::conjunction_v<std::is_convertible<Args, std::string> ...>, void>::type
    addValuesFromCache(Args &&...args) {
        ((data_[args].emplace_back(cache[args]), cache[args] = cache_type{}), ...);
    }

    const std::vector<std::string> &getKeys() {
        return keys_;
    }
    const void clearData() {
        data_ = std::map<std::string, std::vector<cache_type>> {};
    }
    friend std::ostream &operator<<(std::ostream &out, HistogramMeasurement &measurement);

    std::map<std::string, cache_type> cache{};
    static constexpr auto delimeter = ';';
private:
    std::string owner_id_{};
    int number_of_elements_{};
    std::vector<std::string> keys_{};
    std::map<std::string, std::vector<cache_type>> data_;
};
std::ostream &operator<<(std::ostream &out, HistogramMeasurement &measurement);
#endif //CORE_ANALYSER_HISTOGRAM_MEASURMENT_H_
