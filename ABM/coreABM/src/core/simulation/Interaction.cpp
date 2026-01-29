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

#include <sstream>

#include "core/simulation/Interaction.h"
#include "core/analyser/Analyser.h"
#include "core/simulation/states/InteractionState.h"
#include "core/simulation/Cell.h"
#include "core/simulation/factories/InteractionStateFactory.h"
#include "core/simulation/factories/InteractionFactory.h"


Interaction::Interaction(std::string identifier, Cell *cell1, Cell *cell2, double time_delta, double current_time)
        : interactionId(cell1->getSite()->getInteractionFactory()->generateInteractionId()) {
    cellOne = cell1;
    cellTwo = cell2;
    identifier_ = identifier;
    isActiven = true;
    currentCondition = 0;
    setDelete = false;
}

void Interaction::setInitialState(double time_delta, double current_time, Cell *initiatingCell) {

    interactionState = cellOne->getSite()->getInteractionStateFactory()->createInteractionState(this, "InitialInteractionState", cellOne,cellTwo);
    if (cellularConditions.find(initiatingCell) != cellularConditions.end()) {
        currentCondition = cellularConditions[initiatingCell].get();
    }
    interactionState->stateTransition(time_delta, current_time);
}

void Interaction::setState(std::string nameOfState) {
    this->oldinteractionState = this->interactionState;
    this->interactionState = cellOne->getSite()->getInteractionStateFactory()->createInteractionState(this, nameOfState, cellOne, cellTwo);
}

void Interaction::handle(Cell *cell, double timestep, double current_time) {
    interactionState->handleInteraction(cell, timestep, current_time);
}

Cell *Interaction::getFirstCell() {
    return cellOne;
}

Cell *Interaction::getSecondCell() {
    return cellTwo;
}

Cell *Interaction::getOtherCell(Cell *cell) {
    Cell *otherCell;
    if (cell == cellOne) {
        otherCell = cellTwo;
    } else {
        otherCell = cellOne;
    }
    return otherCell;
}

bool Interaction::isActive() {
    return isActiven;
}

void Interaction::close() {
    isActiven = false;
}

void Interaction::fireInteractionEvent(InteractionEvent *ievent, double current_time) {
    cellOne->recieveInteractionEvent(ievent, current_time);
    cellTwo->recieveInteractionEvent(ievent, current_time);
}

std::string Interaction::getInteractionName() const {
    return "Interaction";
}

Condition *Interaction::getCurrentCondition() {
    return currentCondition;
}

std::shared_ptr<Collision> Interaction::getNextCollision() {
    std::shared_ptr<Collision> coll = nullptr;
    if (!currentCollisions.empty()) {
        coll = currentCollisions.front();
        currentCollisions.pop();
    }
    return coll;
}