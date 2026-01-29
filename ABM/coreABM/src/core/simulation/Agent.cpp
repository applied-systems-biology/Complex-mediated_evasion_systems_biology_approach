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

#include <cmath>

#include "core/simulation/Agent.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/states/InteractionState.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"
#include "core/simulation/AgentManager.h"

Agent::Agent() {
    id = 0;
    initialTime = 0;
    is_deleted_ = false;
    initialPosition = std::make_unique<Coordinate3D>();
    previousPosition = std::make_unique<Coordinate3D>();
    movement = std::make_unique<Movement>(site->getNumberOfSpatialDimensions());
    passiveMovement = std::make_unique<Movement>(site->getNumberOfSpatialDimensions());
    hasBeenMoved = false;
    setMovement(movement);
    setPassiveMovement(passiveMovement);
    positionShiftAllowed = true;
    currShift = std::make_unique<Coordinate3D>(Coordinate3D{0, 0, 0});
    setInitialPosition(getPosition());
}

Agent::Agent(std::unique_ptr<Coordinate3D> c, int id, Site *site) {

    this->id = id;
    is_deleted_ = false;
    this->site = site;
    initialPosition = std::make_unique<Coordinate3D>(Coordinate3D{c->x, c->y, c->z});
    previousPosition = std::make_unique<Coordinate3D>();
    setPreviousPosition(initialPosition.get());
    position = std::move(c);
    this->initialTime = 0;
    currShift = 0;
    hasBeenMoved = false;
    positionShiftAllowed = true;
    currShift = std::make_unique<Coordinate3D>(Coordinate3D{0, 0, 0});
    setInitialPosition(getPosition());

}

void Agent::doAllActionsForTimestep(double timestep, double current_time) {
    move(timestep, current_time);
}

void Agent::move(double timestep, double current_time) {
}

int Agent::getId() const {
    return id;
}

Coordinate3D Agent::getPosition() {
    return *position;
}

void Agent::setPosition(Coordinate3D newPos) {
    setPreviousPosition(position.get());
    *position = newPos;
    SphereRepresentation *sphereRep = getMorphology()->getAllSpheresOfThis().front();
    site->getNeighbourhoodLocator()->updateDataStructures(sphereRep);
}

bool Agent::shiftPosition(Coordinate3D *shifter, double current_time, SphereRepresentation *sphereRep, std::string origin) {

    if (positionShiftAllowed) {
        *currShift = *shifter;
        setPreviousPosition(position.get());

        if (sphereRep != 0) {
            sphereRep->shiftPosition(shifter);
        } else {
            sphereRep = getMorphology()->getAllSpheresOfThis().front();
            *position += *shifter;
        }
        if (!site->containsPosition(getPosition())) {
            site->handleBoundaryCross(this, currShift.get(), current_time);
            if (!is_deleted_) {
                site->getNeighbourhoodLocator()->updateDataStructures(sphereRep);
            }
        } else {
            site->getNeighbourhoodLocator()->updateDataStructures(sphereRep);
        }
    }

    return positionShiftAllowed;
}

double Agent::getLifetime(double curTime) {

    return curTime - initialTime;
}

double Agent::getInitialTime() {
    return initialTime;
}

Coordinate3D *Agent::getCurrentShift() {
    return currShift.get();
}

std::string Agent::generatePovObject() {
    return "";
}

Site *Agent::getSite() {
    return site;
}

bool Agent::hasBeenMovedThisTimestep() {
    return hasBeenMoved;
}

void Agent::setBeenMovedThisTimestep(bool newHasBeenMoved) {
    hasBeenMoved = newHasBeenMoved;
    return;
}

void Agent::setInitialPosition(Coordinate3D initPos) {
    *initialPosition = initPos;
}

void Agent::setDeleted() {
    is_deleted_ = true;
}

void Agent::setFeatureValueByName(std::string featureName, int value) {
}

void Agent::applyMethodByName(std::string mehtodName) {
}

void Agent::resetAgent(Coordinate3D pos, double current_time) {
    initialTime = current_time;
    setInitialPosition(pos);
    setPreviousPosition(&pos);
    id = site->getAgentManager()->getIdHandling();
    site->getAgentManager()->incrementIdHandling();
}

bool Agent::coordinateIsInsideAgent(Coordinate3D *coordinate, Agent *agentToInsert) {
    bool isInside = false;
    Coordinate3D center = this->getMorphology()->getBasicSphereOfThis()->getPosition();
    double radius = this->getMorphology()->getBasicSphereOfThis()->getRadius();
    double radiusOfAgentToInsert = agentToInsert->getMorphology()->getBasicSphereOfThis()->getRadius();
    double distanceToCenter = sqrt(
            pow(coordinate->x - center.x, 2) + pow(coordinate->y - center.y, 2) + pow(coordinate->z - center.z, 2));
    if (distanceToCenter < radius - radiusOfAgentToInsert) {
        isInside = true;
    }
    return isInside;
}

Coordinate3D Agent::getCoordinateWithinAgent(Agent *agentToInsert) {
    Coordinate3D pos;
    Coordinate3D center = this->getMorphology()->getBasicSphereOfThis()->getPosition();
    double radius = this->getMorphology()->getBasicSphereOfThis()->getRadius();
    double radiusOfAgentToInsert = agentToInsert->getMorphology()->getBasicSphereOfThis()->getRadius();

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    int rejections = 0;
    do {
        x = site->getRandomGenerator()->generateDouble(center.x - radius + radiusOfAgentToInsert,
                                                       center.x + radius - radiusOfAgentToInsert);
        y = site->getRandomGenerator()->generateDouble(center.y - radius + radiusOfAgentToInsert,
                                                       center.y + radius - radiusOfAgentToInsert);
        z = site->getRandomGenerator()->generateDouble(center.z - radius + radiusOfAgentToInsert,
                                                       center.z + radius - radiusOfAgentToInsert);
        pos = Coordinate3D{x, y, z};
        ++rejections;
    } while (!coordinateIsInsideAgent(&pos, agentToInsert) || !site->containsPosition(pos));

    return pos;
}

bool Agent::hasCollisionsInsideAgent(Agent *ingestedAgent, const std::vector<Agent *> &complex) {
    bool hasCollision = false;
    Coordinate3D pos = ingestedAgent->getMorphology()->getBasicSphereOfThis()->getPosition();
    for (const auto &part : complex) {
        Coordinate3D pos2 = part->getPosition();
        if (part != this && part != ingestedAgent && coordinateIsInsideAgent(&pos2, part)) {
            double distance =
                    pos.calculateEuclidianDistance(
                            part->getMorphology()->getBasicSphereOfThis()->getPosition());
            if (distance < (ingestedAgent->getMorphology()->getBasicSphereOfThis()->getRadius()
                            + part->getMorphology()->getBasicSphereOfThis()->getRadius())) {
                hasCollision = true;
                break;
            }
        }
    }
    return hasCollision;
}

void Agent::setTimestepLastTreatment(double current_time) {
    timestepLastTreatment = current_time;
}

bool Agent::agentTreatedInCurrentTimestep(double current_time) const {
    bool treatedInCurrentTimestep = false;
    if (timestepLastTreatment == current_time) {
        treatedInCurrentTimestep = true;
    }
    return treatedInCurrentTimestep;
}

std::vector<Agent*> Agent::getComplex(){
    const auto interactions = getInteractions()->getAllInteractions();
    std::vector<Agent*> agentComplex;
    if(abm::util::isSubstring("ImmuneCell", this->getTypeName())){
        agentComplex.push_back(this);
        for(size_t i = 0; i < interactions.size(); i++){
            if(interactions[i]->getInteractionName().compare("PhagocyteFungusInteraction") == 0){
                if(interactions[i]->getCurrentState()->getStateName().compare("Phagocytose") == 0 || interactions[i]->getCurrentState()->getStateName().compare("Lysis") == 0){
                    if(this->getTypeName().compare(interactions[i]->getFirstCell()->getTypeName()) == 0){
                        agentComplex.push_back(interactions[i]->getSecondCell());
                    }else{
                        agentComplex.push_back(interactions[i]->getFirstCell());
                    }
                }
            }
        }
    }else if(this->getTypeName().compare("CandidaAlbicans") == 0){
        for(size_t i = 0; i < interactions.size(); i++){
            if(interactions[i]->getInteractionName().compare("PhagocyteFungusInteraction") == 0){
                if(interactions[i]->getCurrentState()->getStateName().compare("Phagocytose") == 0 || interactions[i]->getCurrentState()->getStateName().compare("Lysis") == 0){
                    if(interactions[i]->getSecondCell()->getTypeName().compare("Neutrophil") == 0 || interactions[i]->getSecondCell()->getTypeName().compare("Monocyte") == 0){
                        agentComplex = interactions[i]->getSecondCell()->getComplex();
                    }
                    if(interactions[i]->getFirstCell()->getTypeName().compare("Neutrophil") == 0 || interactions[i]->getFirstCell()->getTypeName().compare("Monocyte") == 0){
                        agentComplex = interactions[i]->getFirstCell()->getComplex();
                    }
                }
            }
        }
    }
    return agentComplex;
}
