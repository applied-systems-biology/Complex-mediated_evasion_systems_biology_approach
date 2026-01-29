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

#ifndef COREABM_UPTAKE_H
#define COREABM_UPTAKE_H

#include "core/simulation/interactiontypes/InteractionType.h"

class Uptake :public InteractionType {
  public:
    // Class for 'Complex Formation' interaction

    Uptake() = default;

    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;

    [[nodiscard]] std::string getTypeName() const final;

  private:
};

#endif // COREABM_UPTAKE_H
