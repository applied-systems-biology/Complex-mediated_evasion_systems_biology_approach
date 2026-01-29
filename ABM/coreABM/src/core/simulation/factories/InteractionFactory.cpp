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

#include "InteractionFactory.h"
#include "core/analyser/Analyser.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interactions.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/cells/interaction/AvoidanceInteraction.h"
#include "core/simulation/cells/interaction/IdenticalCellsInteraction.h"
#include "core/simulation/cells/interaction/NoInteraction.h"
#include "core/simulation/cells/interaction/PhagocyteFungusInteraction.h"

#include "core/utils/macros.h"

InteractionFactory::InteractionFactory(
        const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters, bool use_interactions) {

    if (use_interactions) {
        for (const auto &interaction: interaction_parameters) {
            if (!interaction->cell_conditions.empty()) {
                const auto &cell1 = interaction->cell_conditions[0].first;
                const auto &cell2 = interaction->cell_conditions[1].first;
                interaction_pair_types_[std::make_pair(cell1, cell2)] = interaction->name;
                interaction_pair_types_[std::make_pair(cell2, cell1)] = interaction->name;
                interaction_conditions_[interaction->name] = interaction->cell_conditions;
            }
            interaction_types_[interaction->name] = interaction->type;
        }
    }
}

std::shared_ptr<Interaction> InteractionFactory::createInteraction(double time_delta,
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
        }
        if (interaction != nullptr) {
            interaction->addCurrentCollision(collision);
        }
    }
    return interaction;
}

std::shared_ptr<Interaction> InteractionFactory::createAvoidanceInteraction(Cell *cell_1,
                                                                            std::shared_ptr<Collision> collision,
                                                                            double time_delta,
                                                                            double current_time) {
    std::shared_ptr<Interaction> interaction = nullptr;
    const auto &cell_2 = collision->getCollisionCell();
    if (!(cell_2->isDeleted())) {
        interaction = std::make_shared<AvoidanceInteraction>("AvoidanceInteraction",
                                                             cell_1,
                                                             cell_2,
                                                             time_delta,
                                                             current_time);
        interaction->addCurrentCollision(collision);
    }

    return interaction;
}

std::tuple<std::string, std::string> InteractionFactory::retrieveInteractionIdentifier(Cell *cell_1, Cell *cell_2) {
    std::string interaction_identifier{};
    std::string interaction_type{};

    std::string type_cell_1 = cell_1->getTypeName();
    std::string type_cell_2 = cell_2->getTypeName();
    if (const auto &pair = interaction_pair_types_.find(std::make_pair(type_cell_1, type_cell_2)); pair !=
                                                                                                   interaction_pair_types_.end()) {
        bool condition_cell_1 = true;
        bool condition_cell_2 = true;
        for (const auto &condition : interaction_conditions_[pair->second]) {
            if (type_cell_1 == condition.first) {
                if (!condition.second.empty()) {
                    condition_cell_1 = false;
                    for (const auto &state: condition.second) {
                        if (cell_1->getCurrentCellState()->getStateName() == state) {
                            condition_cell_1 = true;
                        }
                    }
                }
            } else if (type_cell_2 == condition.first) {
                if (!condition.second.empty()) {
                    condition_cell_2 = false;
                    for (const auto &state: condition.second) {
                        if (cell_2->getCurrentCellState()->getStateName() == state) {
                            condition_cell_2 = true;
                        }
                    }
                }
            }
        }
        if (condition_cell_1 && condition_cell_2) {
            interaction_identifier = pair->second;
            interaction_type = interaction_types_[pair->second];
        }
    }
    if (interaction_identifier.empty()) {
        interaction_identifier = "IdenticalCellsInteraction";
        interaction_type = "IdenticalCellsInteraction";
    }
    return std::forward_as_tuple(interaction_type, interaction_identifier);
}

unsigned int InteractionFactory::generateInteractionId() {
    return interaction_id_++;
}

bool InteractionFactory::isInteractionsOn() {
    return !(InteractionFactory::interaction_types_.empty() && InteractionFactory::interaction_pair_types_.empty());
}
