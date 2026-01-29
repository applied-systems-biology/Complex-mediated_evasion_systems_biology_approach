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

#ifndef CORE_SIMULATION_CELLFACTORY_H
#define CORE_SIMULATION_CELLFACTORY_H

#include <string>
#include <array>
#include <map>

#include "core/simulation/Cell.h"
#include "core/utils/io_util.h"


class CellFactory {
public:
  // Factory class for initializing new agents according to the simulator configuration.
    CellFactory(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters);

    virtual std::shared_ptr<Cell> createCell(const std::string &agenttype,
                                            std::unique_ptr<Coordinate3D> c,
                                            int id,
                                            Site *site,
                                            double time_delta,
                                            double current_time, std::string );
    std::map<std::string, std::shared_ptr<abm::util::SimulationParameters::AgentParameters>>
    getAgentConfiguration() { return agent_configurations_;}

protected:
    std::map<std::string, std::shared_ptr<abm::util::SimulationParameters::AgentParameters>> agent_configurations_;
};

#endif /* CORE_SIMULATION_CELLFACTORY_H */