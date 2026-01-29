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

//
// Created by ybachelot on 16.03.23.
//

#include "ComplexFormation.h"
#include "apps/CME/cells/Complex.h"
#include "apps/CME/cells/Drug.h"
#include "apps/CME/interaction/ComplexFormationInteraction.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"

void ComplexFormation::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {
    auto pos = cell->getPosition();
    auto site = cell->getSite();
    auto agent_manager = site->getAgentManager();
    const auto inter = dynamic_cast<ComplexFormationInteraction*>(interaction);
    auto nb_block = 0;
    if (cell->getTypeName() == "AMP"){
        nb_block  = 1;
    }
    else if (cell->getTypeName() == "Drug"){
        auto drug = dynamic_cast<Drug*>(cell);
        nb_block = drug->get_nb_receptors_blocked();
    }
    auto new_comp = agent_manager->emplace_back(site->getCellFactory()->createCell("Complex",
                                                         std::make_unique<Coordinate3D>(pos),
                                                         agent_manager->generateNewID(), site, timestep, current_time, cell->getTypeName()));
    auto comp = std::dynamic_pointer_cast<Complex>(new_comp);
    comp->decrease_free_receptor(nb_block);
    comp->update_radius_diffusion();
    site->getMeasurments()->addValues<HistogramMeasurement>("complex_formation", std::make_pair("time", current_time), std::make_pair("molecule", cell->getTypeName()));
}

std::string ComplexFormation::getTypeName() const {
    return "ComplexFormation";
}