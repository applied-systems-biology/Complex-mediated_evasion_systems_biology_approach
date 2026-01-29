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

#include "IdenticalCellsInteraction.h"


IdenticalCellsInteraction::IdenticalCellsInteraction(std::string identifier,
                                                     Cell *cell1,
                                                     Cell *cell2,
                                                     double time_delta,
                                                     double current_time) : Interaction(identifier,
                                                                                        cell1,
                                                                                        cell2,
                                                                                        time_delta,
                                                                                        current_time) {
    setInitialState(time_delta, current_time, cell1);
}

std::string IdenticalCellsInteraction::getInteractionName() const {
    return "IdenticalCellsInteraction";
}