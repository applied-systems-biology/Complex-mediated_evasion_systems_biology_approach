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

#ifndef CORE_SIMULATION_CONTACTING_H
#define CORE_SIMULATION_CONTACTING_H

#include "InteractionType.h"

class Contacting : public InteractionType {
public:
  // Class for 'Contacting' interaction that is the default interaction between similar cells on collision
    Contacting(bool adhere, double mustOverhead) : adhere(adhere), mustOverhead(mustOverhead) {}

    void handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) final;
    [[nodiscard]] std::string getTypeName() const final;

private:
    bool adhere;
    double mustOverhead;
};

#endif /* CORE_SIMULATION_CONTACTING_H */
