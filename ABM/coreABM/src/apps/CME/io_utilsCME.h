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

#ifndef CORE_UTILS_IO_UTILSCME_H
#define CORE_UTILS_IO_UTILSCME_H

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
#include "apps/CME/Rate/ConcentrationDependentRate.h"

namespace abm::utilCME {

    // To create your own parameters, you have to define them here and let them inherit from the base parameter structs
    // In the Site (here: CuboidSiteCME) you must read in abm::utilCME::getSim.. instead of abm::util:getSim..
    // The parameters are set in the in io_util_CME.cpp
    struct StoppingCriteria{
        std::string molecule;
        double threshold;
        double time;
        double init_time;
    };
    struct CuboidSiteCMEParameters: util::SimulationParameters::CuboidSiteParameters{
        StoppingCriteria stopping_criteria{};
        std::map<std::string, double> flow{};
    };


    struct secretion_rate{
        std::string type;
        std::string mol_uptaken;
        double rate;
        double lag;
        bool spatial;
    };

    struct AgentsParametersCME: util::SimulationParameters::AgentParameters{
        std::vector<std::pair<std::string, secretion_rate>> secretion_;
        bool death_;
    };
    struct Complex: AgentsParametersCME{
        int nb_receptors;
    };
    struct Drug: AgentsParametersCME{
        int nb_recep_blocked;
    };
    struct ConcentrationDependentRateParameters: util::InputParameters::DefaultRateParameters{
        double alpha{};
        ConcentrationDependentRateParameters(const DefaultRateParameters& base)
        : DefaultRateParameters(base){}
    };


    util::SimulationParameters getSimulationParameters(const std::string &simulator_config);
    util::InputParameters getInputParameters(const std::string &input_config);



}

#endif //CORE_UTILS_IO_UTILSCME_H
