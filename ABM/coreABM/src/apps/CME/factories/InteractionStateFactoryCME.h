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

#ifndef INTERACTIONSTATEFACTORYCME_H
#define    INTERACTIONSTATEFACTORYCME_H

#include <variant>

#include "apps/CME/interactiontypes/ComplexExtension.h"
#include "apps/CME/interactiontypes/ComplexFormation.h"
#include "apps/CME/interactiontypes/Uptake.h"
#include "core/simulation/interactiontypes/Contacting.h"
#include "core/simulation/interactiontypes/Ingestion.h"
#include "core/simulation/interactiontypes/RigidContacting.h"
#include "core/simulation/states/InteractionState.h"

class InteractionStateFactoryCME : public InteractionStateFactory {

    using StateSetup = std::map<std::string, std::pair<std::variant<InteractionType, Contacting, RigidContacting, Ingestion, ComplexFormation, ComplexExtension, Uptake>, std::map<std::string, Rate *>>>;

public:
  // Factory class for adding some interactions such as ComplexFormation and Uptake.
  InteractionStateFactoryCME(const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters, RateFactory *rate_factory);

  std::unique_ptr<InteractionState> createInteractionState(Interaction *interaction,
                                                                    std::string interactionStateType,
                                                                    Cell *cell1, Cell *cell2);
protected:

  std::map<std::string, StateSetup> state_parameters_;
};

#endif    /* INTERACTIONSTATEFACTORYCME_H */

