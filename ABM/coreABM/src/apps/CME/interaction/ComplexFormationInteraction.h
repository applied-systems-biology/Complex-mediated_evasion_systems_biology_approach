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

#ifndef COREABM_COMPLEXFORMATIONINTERACTION_H
#define COREABM_COMPLEXFORMATIONINTERACTION_H

#include <string>

#include "core/simulation/Interaction.h"


class ComplexFormationInteraction : public Interaction {

  public:
    // Class for complex formation interaction. This class provides the main functionality if a complex formation event is triggered in a event chain.
    ComplexFormationInteraction(std::string identifier,
                               Cell *cellOne,
                               Cell *cellTwo,
                               double time_delta,
                               double current_time);

    [[nodiscard]] std::string getInteractionName() const final;
};
#endif // COREABM_COMPLEXFORMATIONINTERACTION_H
