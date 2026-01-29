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

#include "ComplexExtension.h"
#include "apps/CME/cells/Amp.h"
#include "apps/CME/interaction/ComplexExtensionInteraction.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"
#include "apps/CME/cells/Complex.h"

void ComplexExtension::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {
    //dynamic_cast<Complex*>(cell)->increase_AMP();
    //SYSTEM_STDOUT("Complex extended!");
    //std::cout << "Now complex with " << dynamic_cast<Complex*>(cell)->get_number_AMP() << " AMPs" << std::endl;
}

std::string ComplexExtension::getTypeName() const {
    return "ComplexExtension";
}