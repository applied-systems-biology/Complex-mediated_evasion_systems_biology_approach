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

#ifndef CORE_UTILS_MACRO_H
#define CORE_UTILS_MACRO_H

#include <iostream>
#include <string>

namespace abm::macro {
    constexpr auto reformulatePrettyFunction(std::string_view &&pretty_function) {
        const auto end = pretty_function.rfind('(');
        const auto begin = pretty_function.substr(0, end).rfind(' ') + 1;
        return pretty_function.substr(begin, end - begin);
    }
}
#define GET_METHOD_NAME_ abm::macro::reformulatePrettyFunction(__PRETTY_FUNCTION__)

#if LOG_ERROR || LOG_INFO || LOG_DEBUG
#define ERROR_STDERR(x) do {constexpr auto _tmp_macro_ = GET_METHOD_NAME_;std::clog << "[ERROR]\t\t"<< '[' << _tmp_macro_ << "] " << x<< '\n';}while (0)
#else
#define ERROR_STDERR(x) do {} while (0)
#endif

#if LOG_INFO || LOG_DEBUG
#define INFO_STDOUT(x) do {constexpr auto _tmp_macro_ = GET_METHOD_NAME_;std::cout << "[INFO]\t\t" << '[' << _tmp_macro_ << "] " << x  << '\n';}while (0)
#else
#define INFO_STDOUT(x) do {} while (0)
#endif

#if LOG_SYSTEM || LOG_DEBUG
#define SYSTEM_STDOUT(x) do {constexpr auto _tmp_macro_ = GET_METHOD_NAME_;std::cout << "[SYSTEM]\t"<< '[' << _tmp_macro_ << "] " << x  << '\n';}while (0)
#else
#define SYSTEM_STDOUT(x) do {} while (0)
#endif

#if LOG_DEBUG
#define DEBUG_STDOUT(x) do {constexpr auto _tmp_macro_ = GET_METHOD_NAME_;std::cout << "[DEBUG]\t\t"<< '[' << _tmp_macro_ << "] " << x  << '\n';}while (0)
#else
#define DEBUG_STDOUT(x) do {} while (0)
#endif


#endif /* CORE_UTILS_MACRO_H */