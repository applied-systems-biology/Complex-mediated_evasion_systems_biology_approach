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

#ifndef CORE_UTILS_IO_UTIL_Example_H
#define CORE_UTILS_IO_UTIL_Example_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <array>
#include <type_traits>

#include <boost/filesystem.hpp>

#include "core/basic/Coordinate3D.h"
#include "core/utils/io_util.h"

namespace abm::utilExample {

    // To create your own parameters, you have to define them here and let them inherit from the base parameter structs
    // In the Site (here: CuboidSiteExample) you must read in abm::utilExample::getSim.. instead of abm::util:getSim..
    // The parameter are set in the in io_util_example.cpp
    struct FungalParametersExample : abm::util::SimulationParameters::FungalParameters {
        bool hyphal_growth{};
        double other_example_parameter{};
    };

    abm::util::SimulationParameters getSimulationParameters(const std::string &simulator_config);

}

#endif //CORE_UTILS_IO_UTILS_Example_H
