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

#include <variant>
#include "core/utils/misc_util.h"

#include "InteractionStateFactory.h"
#include "core/simulation/interactiontypes/Ingestion.h"
#include "core/simulation/interactiontypes/Contacting.h"
#include "core/simulation/interactiontypes/RigidContacting.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/rates/Rate.h"
#include "core/utils/macros.h"
#include "RateFactory.h"

InteractionStateFactory::InteractionStateFactory(const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters,
                                                 RateFactory *rate_factory) {
    for (const auto &parameters: interaction_parameters) {
        for (const auto &state: parameters->states) {
            std::map<std::string, Rate *> state_setup;
            for (const auto&[next_state, rate_name]:state.next_states) {
                state_setup.emplace(next_state, rate_factory->getRate(rate_name));
            }
            if (state.interaction_type == "InteractionType") {
                state_parameters_[parameters->name].emplace(state.name, std::make_pair(InteractionType(),
                                                                                       std::move(state_setup)));
            } else if (state.interaction_type == "Contacting") {
                state_parameters_[parameters->name].emplace(state.name,
                                                            std::make_pair(
                                                                    Contacting(state.adhere, state.must_overhead),
                                                                    std::move(state_setup)));
            } else if (state.interaction_type == "RigidContacting") {
                state_parameters_[parameters->name].emplace(state.name,
                                                            std::make_pair(RigidContacting(state.must_overhead),
                                                                           std::move(state_setup)));
            } else if (state.interaction_type == "Ingestion") {
                state_parameters_[parameters->name].emplace(state.name,
                                                            std::make_pair(Ingestion(), std::move(state_setup)));
            }
        }
    }

}

std::unique_ptr<InteractionState> InteractionStateFactory::createInteractionState(Interaction *interaction,
                                                                                  std::string interactionStateType,
                                                                                  Cell *cell1,
                                                                                  Cell *cell2) {
    std::string identifier = interaction->getIdentifier();
    auto&[type, next_states] = state_parameters_.at(identifier).at(interactionStateType);
    const auto clone = abm::util::overloaded{
            [&](InteractionType type) -> std::unique_ptr<InteractionType> {
                return std::make_unique<InteractionType>(type);
            },
            [&](Contacting type) -> std::unique_ptr<InteractionType> { return std::make_unique<Contacting>(type); },
            [&](RigidContacting type) -> std::unique_ptr<InteractionType> {
                return std::make_unique<RigidContacting>(type);
            },
            [&](Ingestion type) -> std::unique_ptr<InteractionType> { return std::make_unique<Ingestion>(type); },
    };

    std::unique_ptr<InteractionState>
            intState = std::make_unique<InteractionState>(interactionStateType, interaction, std::visit(clone, type),
                                                          next_states.empty());

    if (intState->getInteractionType() == "Ingestion") {
        intState = std::make_unique<InteractionState>(interactionStateType, interaction, std::visit(clone, type),
                                                      false);
        if (abm::util::isSubstring("ImmuneCell", interaction->getFirstCell()->getTypeName())) {
            interaction->getFirstCell()->addIngestions(interaction->getSecondCell()->getId());
        } else if (abm::util::isSubstring("ImmuneCell", interaction->getSecondCell()->getTypeName())) {
            interaction->getSecondCell()->addIngestions(interaction->getFirstCell()->getId());
        }
    }

    intState->addNextStateWithRate(next_states);

    return intState;
}

