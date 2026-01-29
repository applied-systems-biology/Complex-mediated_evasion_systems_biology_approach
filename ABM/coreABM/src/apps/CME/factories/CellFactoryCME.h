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

#ifndef CORE_SIMULATION_CellFactoryCME_H
#define CORE_SIMULATION_CellFactoryCME_H

#include <string>
#include <array>
#include <map>

#include "core/simulation/Cell.h"
#include "core/utils/io_util.h"

class CellFactoryCME : public CellFactory {
public:
  // Factory class for initializing new agents according to the simulator configuration.
    CellFactoryCME(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters> &site_parameters);

    std::shared_ptr<Cell> createCell(const std::string &agenttype,
                                            std::unique_ptr<Coordinate3D> c,
                                            int id,
                                            Site *site,
                                            double time_delta,
                                            double current_time,
                                            std::string mol_type);
};

#endif /* CORE_SIMULATION_CellFactoryCME_h */