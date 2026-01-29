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
// Created by ybachelot on 17.03.23.
//

#include "Uptake.h"

#include "core/simulation/Cell.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"
#include "apps/CME/cells/Pathogenic_cell.h"

void Uptake::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {
    if (cell->getTypeName()=="Pathogenic_cell"){
        dynamic_cast<Pathogenic_cell*>(cell)->update_total_uptake();
        dynamic_cast<Pathogenic_cell*>(cell)->update_current_uptake(current_time, timestep, interaction->getOtherCell(cell)->getTypeName(), interaction->getOtherCell(cell)->getCurrentPosition());
        interaction->getOtherCell(cell)->setDeleted();
        cell->getSite()->getMeasurments()->addValues<HistogramMeasurement>("bounds_distrib", std::make_pair("agent", 1), std::make_pair("nb_bounds", interaction->getOtherCell(cell)->get_nb_bound()));
        cell->getSite()->getMeasurments()->addValues<HistogramMeasurement>("bounds_uptake", std::make_pair("nb_bounds", interaction->getOtherCell(cell)->get_nb_bound()));
    }
    else{
        dynamic_cast<Pathogenic_cell*>(interaction->getOtherCell(cell))->update_total_uptake();
        dynamic_cast<Pathogenic_cell*>(interaction->getOtherCell(cell))->update_current_uptake(current_time, timestep, cell->getTypeName(), cell->getPosition());
        cell->setDeleted();
        cell->getSite()->getMeasurments()->addValues<HistogramMeasurement>("bounds_uptake", std::make_pair("nb_bounds", cell->get_nb_bound()));
        cell->getSite()->getMeasurments()->addValues<HistogramMeasurement>("bounds_distrib", std::make_pair("agent", 1), std::make_pair("nb_bounds", cell->get_nb_bound()));
    }
}

std::string Uptake::getTypeName() const {
    return "Uptake";
}