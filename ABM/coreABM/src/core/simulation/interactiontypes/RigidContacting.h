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

#ifndef CORE_SIMULATION_RIGIDCONTACTING_H
#define CORE_SIMULATION_RIGIDCONTACTING_H

#include "InteractionType.h"
#include "core/simulation/Cell.h"

class RigidContacting : public InteractionType {
public:
    RigidContacting(double mustOverhead) : mustOverhead(mustOverhead) {}
    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;
    std::string getTypeName() const final { return "RigidContacting"; };

private:
    double mustOverhead;
};

#endif /* CORE_SIMULATION_RIGIDCONTACTING_H */
