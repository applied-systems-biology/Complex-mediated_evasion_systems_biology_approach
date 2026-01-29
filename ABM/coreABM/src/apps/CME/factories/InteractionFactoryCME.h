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

#ifndef INTERACTIONFACTORYCME_H
#define    INTERACTIONFACTORYCME_H

#include <memory>
#include <utility>
#include <map>

#include "core/simulation/factories/InteractionFactory.h"
#include "core/simulation/Interaction.h"
#include "core/analyser/Analyser.h"

class Collision;
class Cell;

class InteractionFactoryCME : public InteractionFactory {
public:
  // Factory class for the interactions specified in the simulator configuration
    InteractionFactoryCME(const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters, bool use_interactions);

    std::shared_ptr<Interaction> createInteraction(double time_delta,
                                                          double current_time,
                                                          const std::shared_ptr<Collision> &collision,
                                                          InSituMeasurements *measurements) override;

};

#endif    /* INTERACTIONFACTORYCME_H */

