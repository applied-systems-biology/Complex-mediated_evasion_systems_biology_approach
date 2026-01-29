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
// Created by ybachelot on 15.01.25.
//

#include "ConcentrationDependentRate.h"
#include "apps/CME/cells/Pathogenic_cell.h"

double ConcentrationDependentRate::calculateProbability(double timestep, Condition *cond,  Cell* cell, Site* site) {
    Pathogenic_cell *pathogen = dynamic_cast<Pathogenic_cell*>(cell);
    auto lower_limit = site->getLowerLimits();
    auto upper_limit = site->getUpperLimits();
    auto volume = ((upper_limit.x - lower_limit.x) * (upper_limit.y - lower_limit.y) * (upper_limit.z - lower_limit.z)) - cell->getMorphology()->getVolume();

    auto inside_conc = pathogen->get_total_uptake() / volume;
    current_rate_ = initial_rate_ * (1.0 + (alpha_ *  inside_conc));
    return current_rate_ * timestep;
}