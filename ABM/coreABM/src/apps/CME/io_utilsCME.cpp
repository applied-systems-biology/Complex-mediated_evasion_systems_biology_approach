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

#include <map>
#include <optional>
#include <random>

#include "io_utilsCME.h"
#include "external/json.hpp"
#include "core/utils/macros.h"

using json = nlohmann::json;

abm::util::SimulationParameters abm::utilCME::getSimulationParameters(const std::string &simulator_config){
    // Load parameters with base abm::utils
    auto simp = abm::util::getSimulationParameters(simulator_config);

    std::ifstream json_file(simulator_config);
    json json_parameters;
    json_file >> json_parameters;
    // Load CuboidWBIA parameters
    auto& sitep = simp.site_parameters;
    auto cs_para = CuboidSiteCMEParameters{*sitep};
    for (auto& site : json_parameters["Agent-Based-Framework"]["Sites"]) {
        cs_para.lower_bound = {site["CuboidSite"]["x_range"][0], site["CuboidSite"]["y_range"][0], site["CuboidSite"]["z_range"][0]};
        cs_para.upper_bound = {site["CuboidSite"]["x_range"][1], site["CuboidSite"]["y_range"][1], site["CuboidSite"]["z_range"][1]};
        for (const auto &flow: site["CuboidSite"]["flow"]) {
            auto name = flow["agent"];
            auto rate = flow.value("rate", 0);
            cs_para.flow[name] = rate;
        }
        for (const auto &stopping: site["CuboidSite"]["stopping_criteria"]) {
            StoppingCriteria stop_param;
            stop_param.molecule = stopping["agent"];
            stop_param.threshold = stopping.value("threshold", 0.0);
            stop_param.time = stopping.value("time", 1000000);
            stop_param.init_time = stop_param.time;
            cs_para.stopping_criteria = stop_param;
        }
        for (auto& agent : cs_para.agent_manager_parameters.agents) {
            for (const auto& agent_type : site["AgentManager"]["Types"]) {
                if (agent->type == agent_type) {
                    if (agent->type == "Complex") {
                        auto agent_ = site["AgentManager"]["Agents"].at(std::string(agent_type));
                        auto ag_comp = Complex{*agent};
                        ag_comp.nb_receptors = agent_.value("nb_receptor", 1);
                        agent = std::make_shared<Complex>(ag_comp);
                    }
                    else if (agent->type == "Drug"){
                        auto agent_ = site["AgentManager"]["Agents"].at(std::string(agent_type));
                        auto ag_drug = Drug{*agent};
                        ag_drug.nb_recep_blocked = agent_.value("nb_receptor_blocked", 1);
                        agent = std::make_shared<Drug>(ag_drug);
                    }
                    else {
                        auto agent_ = site["AgentManager"]["Agents"].at(std::string(agent_type));
                        auto ag_para = AgentsParametersCME{*agent};
                        ag_para.death_ = agent_.value("death", false);
                        for (const auto& mol : agent_["Secretion"]) {
                            auto secretion = agent_.at(std::string(mol));
                            const auto type_ = secretion.value("type", "constant");
                            std::unique_ptr<secretion_rate> sec_rate = std::make_unique<secretion_rate>();
                            sec_rate->type = type_;
                            sec_rate->rate = secretion.value("rate", 0.0);
                            sec_rate->spatial = secretion.value("spatial", false);
                            sec_rate->lag = secretion.value("lag", 0.0);
                            if (type_ == "variable") {
                                sec_rate->mol_uptaken = secretion.value("mol_uptaken", "");
                            }
                            ag_para.secretion_.emplace_back(std::make_pair(std::string(mol), *sec_rate));
                        }
                        agent = std::make_shared<AgentsParametersCME>(ag_para);
                    }
                }
            }
        }
        sitep = std::make_unique<CuboidSiteCMEParameters>(cs_para);
    }
    json_file.close();
    return std::move(simp);
}

abm::util::InputParameters abm::utilCME::getInputParameters(const std::string &input_config) {
    auto input_param = abm::util::getInputParameters(input_config);

    // Load file again to load CME specific parameters
    std::ifstream json_file(input_config);
    json json_parameters;
    json_file >> json_parameters;
    for (auto& rates: input_param.rates) {
        for (auto& rate : json_parameters["Agent-Based-Framework"]["Rates"].items()) {
            auto rate_type_ = "ConcentrationDependentRate";
            if(rates->type == rate_type_){
                if (rate.value().value("type", "ConstantRate") ==  rate_type_) {
                    auto conc_rate = ConcentrationDependentRateParameters{*rates};
                    //DEBUG_STDOUT(rate.value());
                    conc_rate.alpha = rate.value().value("alpha", 0.0);
                    rates = std::make_unique<ConcentrationDependentRateParameters>(conc_rate);
                }
            }
        }
    }
    json_file.close();
    return std::move(input_param);
}
