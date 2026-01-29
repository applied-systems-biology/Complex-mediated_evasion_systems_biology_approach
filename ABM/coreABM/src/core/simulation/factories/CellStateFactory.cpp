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

#include "CellStateFactory.h"

#include "RateFactory.h"
#include "core/utils/macros.h"
#include "core/simulation/Cell.h"

CellStateFactory::CellStateFactory(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters,
                                   RateFactory *rate_factory) {
    for (const auto &agent : site_parameters->agent_manager_parameters.agents) {
        StateSetup state_setup{};
        for (const auto &state : agent->states) {
            state_setup[state.name] = {};
            for (const auto&[next_state, rate_name] : state.next_states) {
                state_setup[state.name].emplace(next_state, rate_factory->getRate(rate_name));
            }
        }
        next_states[std::make_pair(site_parameters->identifier, agent->type)] = state_setup;
    }
}

std::shared_ptr<CellState> CellStateFactory::createCellState(Cell *cell, const std::string &state_name) {
    return std::make_shared<CellState>(state_name,
                                       cell,
                                       CellStateFactory::getNextStates(state_name, cell->getSite()->getIdentifier(),
                                                                       cell->getTypeName()));
}

std::map<std::string, Rate *> CellStateFactory::getNextStates(const std::string &state_name,
                                                                    const std::string &site_identifier,
                                                                    const std::string &agent_type) {
    if (auto result_pair = next_states.find(std::make_pair(site_identifier, agent_type)); result_pair !=
                                                                                          next_states.end()) {
        if (auto state_pair = result_pair->second.find(state_name); state_pair != result_pair->second.end()) {
            return state_pair->second;
        }
        return {};
    }
    ERROR_STDERR("Agent Type " << agent_type << " does not exist.");
    exit(1);
}

CellStateFactory::StateSetup &
CellStateFactory::getStateSetup(const std::string &site_identifier, const std::string &agent_type) {
    if (auto result_pair = next_states.find(std::make_pair(site_identifier, agent_type)); result_pair !=
                                                                                          next_states.end()) {
        return result_pair->second;
    }
    ERROR_STDERR("Agent Type " << agent_type << " does not exist.");
    exit(1);
}
