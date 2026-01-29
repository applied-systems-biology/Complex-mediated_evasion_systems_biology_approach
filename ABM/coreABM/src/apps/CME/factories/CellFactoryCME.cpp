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

#include "CellFactoryCME.h"
#include "apps/CME/cells/Amp.h"
#include "apps/CME/cells/Complex.h"
#include "apps/CME/cells/Defensive.h"
#include "apps/CME/cells/Drug.h"
#include "apps/CME/cells/Pathogenic_cell.h"
#include "core/simulation/Site.h"
#include "core/simulation/cells/FungalCell.h"
#include "core/simulation/cells/ImmuneCell.h"

CellFactoryCME::CellFactoryCME(const std::unique_ptr<abm::util::SimulationParameters::SiteParameters>
        &site_parameters) : CellFactory(site_parameters){
}

std::shared_ptr<Cell> CellFactoryCME::createCell(const std::string &agenttype,
                                              std::unique_ptr<Coordinate3D> c,
                                              int id,
                                              Site *site,
                                              double time_delta,
                                              double current_time,
                                              std::string mol_type) {
    std::shared_ptr<Cell> agent{};
    auto site_tag = site->getIdentifier();
    if (agenttype == "Cell") {
        agent = std::make_shared<Cell>(std::move(c), id, site, time_delta, current_time);
    } else if (agenttype == "ImmuneCell") {
        agent = std::make_shared<ImmuneCell>(std::move(c), id, site, time_delta, current_time);
    } else if (agenttype == "FungalCell") {
        agent = std::make_shared<FungalCell>(std::move(c), id, site, time_delta, current_time);
    } else if (agenttype == "AMP") {
        agent = std::make_shared<Amp>(std::move(c), id, site, time_delta, current_time);
        //DEBUG_STDOUT(agenttype << " cell created");
    } else if (agenttype == "Drug") {
        agent = std::make_shared<Drug>(std::move(c), id, site, time_delta, current_time);
        //DEBUG_STDOUT(agenttype << " cell created");
    } else if (agenttype == "Complex") {
        agent = std::make_shared<Complex>(std::move(c), id, site, time_delta, current_time);
        const auto ag = std::dynamic_pointer_cast<Complex>(agent);
        ag->add_molecule_to_complex(mol_type);
    } else if (agenttype == "Defensive") {
        agent = std::make_shared<Defensive>(std::move(c), id, site, time_delta, current_time);
        //DEBUG_STDOUT(agenttype << " cell created");
    } else if (agenttype == "Pathogenic_cell") {
        agent = std::make_shared<Pathogenic_cell>(std::move(c), id, site, time_delta, current_time);
    }
    agent->setup(time_delta, current_time, agent_configurations_[site_tag + agenttype].get());
    return agent;
}


