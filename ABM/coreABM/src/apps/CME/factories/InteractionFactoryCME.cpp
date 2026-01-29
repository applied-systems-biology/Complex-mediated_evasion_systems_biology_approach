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

#include "InteractionFactoryCME.h"
#include "core/analyser/Analyser.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/cells/interaction/AvoidanceInteraction.h"
#include "apps/CME/interaction/ComplexFormationInteraction.h"
#include "core/simulation/cells/interaction/IdenticalCellsInteraction.h"
#include "core/simulation/cells/interaction/NoInteraction.h"
#include "core/simulation/cells/interaction/PhagocyteFungusInteraction.h"
#include "apps/CME/interaction/UptakeInteraction.h"

#include "apps/CME/cells/Complex.h"
#include "apps/CME/interaction/ComplexExtensionInteraction.h"
#include "core/utils/macros.h"
#include "apps/CME/cells/Drug.h"

InteractionFactoryCME::InteractionFactoryCME(const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters, bool use_interactions)
    : InteractionFactory(interaction_parameters, use_interactions){}


std::shared_ptr<Interaction> InteractionFactoryCME::createInteraction(double time_delta,
                                                                   double current_time,
                                                                   const std::shared_ptr<Collision> &collision,
                                                                   InSituMeasurements *measurements) {
    std::shared_ptr<Interaction> interaction = nullptr;
    const auto &cell_1 = collision->getCell();
    const auto &cell_2 = collision->getCollisionCell();
    if (!(cell_2->isDeleted())) {
        const auto[interaction_name, identifier] = retrieveInteractionIdentifier(cell_1, cell_2);
        if (interaction_name == "IdenticalCellsInteraction") {
            interaction = std::make_shared<IdenticalCellsInteraction>(identifier, cell_1, cell_2, time_delta,
                                                                      current_time);
        } else if (interaction_name == "NoInteraction") {
            interaction = std::make_shared<NoInteraction>(identifier, cell_1, cell_2, time_delta, current_time);
        } else if (interaction_name == "PhagocyteFungusInteraction") {
            interaction = std::make_shared<PhagocyteFungusInteraction>(identifier, cell_1, cell_2, time_delta,
                                                                       current_time);
        } else if (interaction_name == "ComplexFormationInteraction"){
            if (cell_1->getTypeName() == "AMP" || cell_1->getTypeName() == "Drug") {
                interaction = std::make_shared<ComplexFormationInteraction>(identifier, cell_1, cell_2, time_delta,
                                                                            current_time);
            }
        } else if (interaction_name == "ComplexExtensionInteraction"){
            if (cell_1->getTypeName() == "Complex") {
                auto comp = dynamic_cast<Complex*>(cell_1);
                auto nb_recep_block = 1;
                if (cell_2->getTypeName() == "Drug"){
                    auto drug = dynamic_cast<Drug*>(cell_2);
                    nb_recep_block = drug->get_nb_receptors_blocked();
                }
                if (comp->get_free_receptors() >= nb_recep_block) {
                    interaction = std::make_shared<ComplexExtensionInteraction>(identifier, cell_1, cell_2, time_delta,
                                                                                current_time, nb_recep_block);
                }
            }
        }else if (interaction_name == "UptakeInteraction"){
            if (cell_1->getTypeName() == "AMP") {
                interaction = std::make_shared<UptakeInteraction>(identifier, cell_1, cell_2, time_delta,
                                                                  current_time);
            }
        }
        if (interaction != nullptr) {
            interaction->addCurrentCollision(collision);
        }
    }
    return interaction;
}
