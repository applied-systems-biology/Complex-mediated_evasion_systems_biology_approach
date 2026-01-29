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

#include <optional>
#include <random>

#include "io_utils_example.h"
#include "external/json.hpp"
#include "core/utils/macros.h"

using json = nlohmann::json;

abm::util::SimulationParameters abm::utilExample::getSimulationParameters(const std::string &simulator_config){
    // Load parameters with base abm::utils
    auto simp = abm::util::getSimulationParameters(simulator_config);

    // Load file again to load example specific parameters
    std::ifstream json_file(simulator_config);
    json json_parameters;
    json_file >> json_parameters;

    // write example specific parameters into program
    for (auto &agent: simp.site_parameters->agent_manager_parameters.agents) {
        auto agent_type_ = "FungalCellExample";
        for (const auto &site: json_parameters["Agent-Based-Framework"]["Sites"]) {
            if (agent->type == agent_type_) {
                auto agent_ = site["AgentManager"]["Agents"].at(std::string(agent_type_));
                auto fge = FungalParametersExample{*agent};
                fge.hyphal_growth = agent_.value("hyphal_growth", true);
                fge.other_example_parameter = agent_.value("some_para", 1.0);
                agent = std::make_shared<FungalParametersExample>(fge);
            }

        }
    }
    return std::move(simp);
};