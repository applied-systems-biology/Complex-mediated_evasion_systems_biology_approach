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

#include "core/simulation/factories/InteractionStateFactory.h"
#include "ComplexExtensionInteraction.h"
#include "core/analyser/Analyser.h"
#include "core/simulation/Cell.h"
#include "apps/CME/cells/Complex.h"


ComplexExtensionInteraction::ComplexExtensionInteraction(std::string identifier,
                                                         Cell *cell1,
                                                         Cell *cell2,
                                                         double time_delta,
                                                         double current_time,
                                                         int nb_recep_block) : Interaction(identifier,
                                                                                          cell1,
                                                                                          cell2,
                                                                                          time_delta,
                                                                                          current_time) {
    // Here cell1 is always the complex
    auto complex = dynamic_cast<Complex*>(cell1);
    complex->extension(*cell2, nb_recep_block);
    cell2->setDeleted();
    setInitialState(time_delta, current_time, cell1);
}

std::string ComplexExtensionInteraction::getInteractionName() const {
    return "ComplexExtensionInteraction";
}