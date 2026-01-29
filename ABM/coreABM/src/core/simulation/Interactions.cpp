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

#include <algorithm>

#include "core/simulation/Interactions.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/factories/InteractionFactory.h"
#include "core/simulation/states/InteractionState.h"
#include "core/simulation/neighbourhood/Collision.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/utils/macros.h"


Interactions::Interactions(Cell *cell, NeighbourhoodLocator *nhLocator) {
    this->cell = cell;
    neighbourhoodLocator = nhLocator;
}

Interactions::~Interactions() {
    removeAllInteractions();
}

void Interactions::doWholeProcess(double time_delta, double current_time, InSituMeasurements *measurments) {
    if (cell->getSite()->getInteractionFactory()->isInteractionsOn()) {
        const auto collisions = neighbourhoodLocator->getCollisions(cell);
        for (const auto &collision: collisions) {
            if (auto itCellInteraction = interactionPartners.find(collision->getCollisionCell()); itCellInteraction
                                                                                                  !=
                                                                                                  interactionPartners.end()) {
                collision->type = MeasurementType::EXISTING_INTERACTION;
                itCellInteraction->second->addCurrentCollision(collision);
            } else if (auto interaction = cell->getSite()->getInteractionFactory()->createInteraction(time_delta, current_time, collision,
                                                                                measurments); interaction
                                                                                              != nullptr) {
                collision->type = MeasurementType::NEW_INTERACTION;
                addInteraction(interaction);
                collision->getCollisionCell()->getInteractions()->addInteraction(interaction);
            }
        }
        executeAllInteractions(time_delta, current_time);
    }
}

void Interactions::removeCollisionsOfExistingInteractions(std::vector<std::shared_ptr<Collision>> &collisions) {

    auto itCollisions = collisions.begin();
    while (itCollisions != collisions.end()) {

        auto collision = *itCollisions;
        Cell *collCell = collision->getCollisionCell();
        if (auto itCellInteraction = interactionPartners.find(collCell);itCellInteraction !=
                                                                        interactionPartners.end()) {

            collisions.erase(itCollisions);
        } else {
            itCollisions++;
        }
    }
}

void Interactions::addInteraction(std::shared_ptr<Interaction> interaction) {
    interactions.push_back(interaction);
    Cell *otherCell = interaction->getOtherCell(cell);
    interactionPartners[otherCell] = interaction;
}

void Interactions::executeAllInteractions(double timestep, double current_time) {
    size_t numberOfInteractions = interactions.size();
    if (numberOfInteractions > 0) {
        std::vector<unsigned int> currentPermutation =
                abm::util::generateRandomPermutation(cell->getSite()->getRandomGenerator(), numberOfInteractions);

        for (unsigned int currentInteraction : currentPermutation) {
            auto interaction = interactions.at(currentInteraction);
            if (!(interaction->isDelted())) {
                std::string intName = interaction->getInteractionName();
                interaction->handle(cell, timestep, current_time);
                if (interaction->getOtherCell(cell)->getCurrentCellState()->checkForDeath(current_time)) {
                    continue;
                }
                if (cell->getCurrentCellState()->checkForDeath(current_time) || cell->isDeleted()) {
                    break;
                }
            } else {
                interaction->close();
            }
        }
    }
    removeClosedInteractions();
}

void Interactions::removeClosedInteractions() {
    auto it = interactions.begin();
    while (it != interactions.end()) {
        auto interaction = *it;
        if (!interaction->isActive()) {
            Cell *otherCell = interaction->getOtherCell(cell);
            otherCell->getInteractions()->removeInteraction(interaction.get());
            removeInteraction(interaction.get());
            interaction = nullptr;
        } else {
            it++;
        }
    }
}

void Interactions::removeInteraction(Interaction *interaction) {
    Cell *otherCell = interaction->getOtherCell(cell);

    if (!interactionPartners.empty()) {
        if (interactionPartners.count(otherCell) > 0) {
            interactionPartners.erase(interactionPartners.find(otherCell));
        }
    }
    interactions.erase(std::remove_if(interactions.begin(),
                                      interactions.end(),
                                      [interaction](const auto &i) { return i.get() == interaction; }),
                       interactions.end());
}

void Interactions::removeAllInteractions() {
    auto it = interactions.begin();
    while (it != interactions.end()) {
        auto interaction = *it;
        Cell *otherCell = interaction->getOtherCell(cell);
        otherCell->getInteractions()->removeInteraction(interaction.get());
        removeInteraction(interaction.get());
        interaction = nullptr;
    }
}

void Interactions::displayInteractions() {
    auto it = interactions.begin();
    DEBUG_STDOUT("Interactions of cell-id " + std::to_string(cell->getId()) + ":");
    if (it == interactions.end()) {
        DEBUG_STDOUT(" none");
    }
    while (it != interactions.end()) {
        auto interaction = *it;
        DEBUG_STDOUT(
                interaction->getInteractionName() + " " +
                interaction->getFirstCell()->getCurrentCellState()->getStateName() +
                " " + std::to_string(interaction->getFirstCell()->getId()) + " vs. " +
                interaction->getSecondCell()->getCurrentCellState()->getStateName() +
                " " + std::to_string(interaction->getSecondCell()->getId()));
        it++;
    }
}

bool Interactions::hasInteractions() {
    return !interactions.empty();
}

bool Interactions::hasPhagFungInteraction() {
    bool phagFungInt = false;
    for (const auto &i : interactions) {
        if (i->getInteractionName().compare("PhagocyteFungusInteraction") == 0) {
            phagFungInt = true;
        }
        break;
    }
    return phagFungInt;
}

bool Interactions::hasCollisions() {
    return neighbourhoodLocator->hasCollision(cell);
}

void Interactions::avoidNewInteractions(double time_delta, double current_time) {
    if (cell->getSite()->getInteractionFactory()->isInteractionsOn()) {
        auto collisions = neighbourhoodLocator->getCollisions(cell);
        removeCollisionsOfExistingInteractions(collisions);
        doAvoidanceInteractions(collisions, time_delta, current_time);
    }
}

void Interactions::doAvoidanceInteractions(std::vector<std::shared_ptr<Collision>> &collisions, double time_delta,
                                           double current_time) {
    auto itCollisions = collisions.begin();
    for (const auto &collision: collisions) {

        std::shared_ptr<Interaction> interaction = cell->getSite()->getInteractionFactory()->createAvoidanceInteraction(cell, collision,
                                                                                                  time_delta,
                                                                                                  current_time);
        if (interaction != nullptr) {
            //add interaction procedure and handling
            addInteraction(interaction);
            collision->getCollisionCell()->getInteractions()->addInteraction(interaction);
            interaction->handle(cell, time_delta, current_time);
        }
    }
    removeClosedInteractions();
}

const std::vector<std::shared_ptr<Interaction>> &Interactions::getAllInteractions() {
    return interactions;
}

const std::map<Cell *, std::shared_ptr<Interaction>> &Interactions::getAllInteractionPartners() {
    return interactionPartners;
}


