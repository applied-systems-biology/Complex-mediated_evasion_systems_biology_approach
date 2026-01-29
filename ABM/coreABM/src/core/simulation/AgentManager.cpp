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

#include "apps/CME/cells/Complex.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/AgentManager.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"
#include "core/simulation/boundary-condition/PeriodicBoundaries.h"
#include "core/simulation/cells/FungalCell.h"
#include "core/simulation/cells/interaction/PhagocyteFungusInteraction.h"
#include "core/simulation/factories/CellFactory.h"
#include "core/simulation/factories/RateFactory.h"
#include "core/simulation/states/InteractionState.h"
#include "core/utils/macros.h"

AgentManager::AgentManager(abm::util::SimulationParameters& parameters, Site *site) {
    this->site = site;
    time_delta_ = parameters.time_stepping;
    idHandling = 0;
    idHandlingSphereRepresentation = 0;
    for (auto& agent: parameters.site_parameters->agent_manager_parameters.agents){
        agent_types_.emplace_back(agent->type);
    }
}

void AgentManager::inputOfAgents(double current_time, Randomizer *random_generator) {
    for (auto agent_name: agent_types_) {
        double lambda = site->getInputRate(agent_name);
        if (lambda > 0) {
            if (current_time == 0) {
                // Generate the initial nextInputEventTime
                lastInputEventTime = 0;
                nextInputEventTime = 1.0 / (lambda) * log(1.0 / random_generator->generateDouble());
                DEBUG_STDOUT("Next input of " + agent_name + " at " + std::to_string(nextInputEventTime) + " with lambda input rate =" +
                             std::to_string(lambda));
            }
            if (current_time >= nextInputEventTime) {
                //Generate nextInputEventTime when time the previously drawn nextInputEventTime is exceeded
                insertAgentAtBoundary(site, agent_name, current_time);
                lastInputEventTime = nextInputEventTime;
                double diff = 1.0 / (lambda) * log(1.0 / random_generator->generateDouble());
                nextInputEventTime = lastInputEventTime + diff;
                DEBUG_STDOUT("Insert " + agent_name + " at " + std::to_string(current_time) + ". Next input at "
                + std::to_string(nextInputEventTime) + " with lambda input rate =" + std::to_string(lambda));
            }
        }
    }
}

void AgentManager::insertAgentAtBoundary(Site *site, std::string agentType, double current_time) {
    Coordinate3D initialPosition;
    Coordinate3D initialVector;

    Agent *agent = nullptr;
    int rejections = 0;

    // Add an agent to the system without having any collisions
    do {
        if (agent != nullptr) {
            removeAgent(site, agent, current_time);
        }
        initialPosition = site->getRandomBoundaryPoint();
        initialVector = site->getBoundaryInputVector();
        agent = createAgent(site, agentType, initialPosition, &initialVector, current_time);
        rejections++;
    } while ((agent->getInteractions()->hasCollisions() && rejections < 10000));
    if (rejections > 9999) {
        DEBUG_STDOUT("Could not find a position for the agent at the boundary. Agent is not added to the system.");
        if (agent != nullptr) removeAgent(site, agent, current_time);
    }
}

Agent *AgentManager::createAgent(Site *site, std::string agenttype, Coordinate3D c, Coordinate3D *prevMove, double current_time) {
    auto agent = emplace_back(site->getCellFactory()->createCell(agenttype, std::make_unique<Coordinate3D>(c), generateNewID(),
                                                      site, time_delta_, current_time, ""));
    agent->getMovement()->setPreviousMove(prevMove);

    return agent.get();
}

void AgentManager::replaceAgent(Site *site,Agent *agent, std::unique_ptr<Coordinate3D> newCoord, Coordinate3D *prevMove,
                                double current_time) {

    std::string agentType = agent->getTypeName();
    auto newAgent = site->getCellFactory()->createCell(agentType, std::move(newCoord), idHandling, site, time_delta_, current_time, "");
    if (newAgent != 0) {
        newAgent->getMovement()->setPreviousMove(prevMove);
        std::replace_if(allAgents.begin(), allAgents.end(), [agent](const auto &a) { return agent == a.get(); },
                        newAgent);
        idHandling++;
    }
    if (agentType == "Complex") {
        auto old_ag_comp = dynamic_cast<Complex*>(agent);
        auto new_ag_comp = std::dynamic_pointer_cast<Complex>(newAgent);
        for (auto mol : old_ag_comp->get_complex_content()) {
            new_ag_comp->add_molecule_to_complex(mol);
        }
        new_ag_comp->set_free_receptors(old_ag_comp->get_free_receptors());
    }
    agent->setDeleted();
}

void
AgentManager::replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time) {

    // If the new agent is the nullptr, delete the agent
    if (!newAgent) {
        auto all_interactions = agentToReplace->getInteractions()->getAllInteractions();
        for (size_t i = 0; i < all_interactions.size(); i++) {
            if (all_interactions.at(i)->getInteractionName() == "PhagocyteFungusInteraction") {
                Cell *cell1 = all_interactions.at(i)->getFirstCell();
                Cell *cell2 = all_interactions.at(i)->getSecondCell();
                if (cell1->getCurrentCellState()->getStateName() == "FungalPhagocytosed") {
                    cell1->setDeleted();
                } else if (cell2->getCurrentCellState()->getStateName() == "FungalPhagocytosed"){
                    cell2->setDeleted();
                }
            }
        }
        agentToReplace->setDeleted();
    } else {
        std::replace_if(allAgents.begin(),
                        allAgents.end(),
                        [agent = agentToReplace](const auto &a) { return agent == a.get(); },
                        newAgent);
    }

}

void AgentManager::removeAgent(Site *site, Agent *agent, double current_time) {

    // Remove agent and all its corresponding spheres in the neighbourhoodlocator
    for (auto sphRep: agent->getSurface()->getAllSpheresOfThis()) {
        site->getNeighbourhoodLocator()->removeSphereRepresentation(sphRep);
    }
    agent->setDeleted();
    if (abm::util::isSubstring("FungalCell", agent->getTypeName())) {
        removeFungalCellFromList(agent->getId(), current_time);
    }
    allAgents.erase(std::remove_if(allAgents.begin(),
                                   allAgents.end(),
                                   [agent](const auto &a) { return agent == a.get(); }), allAgents.end());
}

int AgentManager::getAgentQuantity(std::string agenttype) {
    int count = 0;
    for (auto agent: allAgents) {
        if (agent->getCurrentCellState()->getStateName() != "Death" && agent != 0) {
            if (agent->getTypeName() == agenttype) {
                count++;
            }
        }
    }
    return count;
}

const std::vector<std::shared_ptr<Agent>> &AgentManager::getAllAgents() {
    return allAgents;

}

void AgentManager::cleanUpAgents(double current_time) {

    for (auto it = allAgents.begin(); it != allAgents.end();) {
        if (*it == 0) {
            it = allAgents.erase(it);
        } else {
            if ((*it)->isDeleted()) {
                for (const auto &sphere: (*it)->getMorphology()->getAllSpheresOfThis()) {
                    site->getNeighbourhoodLocator()->removeSphereRepresentation(sphere);
                    removeSphereRepresentation(sphere);
                }
                if (abm::util::isSubstring("FungalCell", (*it)->getTypeName())) { this->removeFungalCellFromList((*it)->getId(), current_time); }
                it = allAgents.erase(it);
            } else {
                it++;
            }
        }
    }
}

int AgentManager::getNextSphereRepresentationId(SphereRepresentation *sphereRep) {
    sphereIdToCell[idHandlingSphereRepresentation] =
            sphereRep->getMorphologyElementThisBelongsTo()->getMorphologyThisBelongsTo()->getCellThisBelongsTo();
    sphereIdToSphereRep[idHandlingSphereRepresentation] = sphereRep;
    allSphereRepresentations.insert(sphereRep);
    return idHandlingSphereRepresentation++;
}

Cell *AgentManager::getCellBySphereRepId(int sphereRepId) {
    if (sphereIdToCell.find(sphereRepId) != sphereIdToCell.end()) {
        return sphereIdToCell[sphereRepId];
    }
    return nullptr;
}

SphereRepresentation *AgentManager::getSphereRepBySphereRepId(int sphereRepId) {
    if (sphereIdToSphereRep.find(sphereRepId) != sphereIdToSphereRep.end()) {
        return sphereIdToSphereRep[sphereRepId];
    }
    return nullptr;
}

void AgentManager::removeSphereRepresentation(SphereRepresentation *sphereRep) {
    allSphereRepresentations.erase(sphereRep);
    sphereIdToSphereRep.erase(sphereRep->getId());
    sphereIdToCell.erase(sphereRep->getId());
}

int AgentManager::getIdHandling() const {
    return idHandling;
}

void AgentManager::incrementIdHandling() {
    idHandling++;
}

std::vector<std::string> AgentManager::getAllAgentTypes(){
    std::vector<std::string> agent_names;
    for (auto agent: allAgents) {
        auto name = agent->getTypeName();
        auto it = find(agent_names.begin(), agent_names.end(), name);
        if (it == agent_names.end()){
            agent_names.emplace_back(name);
        }
    }
    return agent_names;
}

double AgentManager::getOccupancyDensityOfSpace() {
    double cellVolume = 0;
    Coordinate3D lowerLim = site->getLowerLimits();
    Coordinate3D upperLim = site->getUpperLimits();

    double x = abs(lowerLim.x) + upperLim.x;
    double y = abs(lowerLim.y) + upperLim.y;
    double z = abs(lowerLim.z) + upperLim.z;
    double siteVolume = 0;
    if (site->getNumberOfSpatialDimensions() == 2) {
        siteVolume = x * y;
    } else if (site->getNumberOfSpatialDimensions() == 3) {
        siteVolume = x * y * z;
    }

    for (size_t i = 0; i < allAgents.size(); i++) {
        cellVolume = cellVolume + allAgents[i]->getMorphology()->getVolume();
    }
    double occupancy = (cellVolume * 100) / siteVolume;
    return occupancy;
}

void AgentManager::removeFungalCellFromList(int id, double current_time) {
    for (size_t i = 0; i < activeFungalCells.size(); i++) {
        if (activeFungalCells.at(i)->getId() == id) {
            activeFungalCells.erase(activeFungalCells.begin() + i);
            fungalCellRemoveTimes.push_back(current_time);
            fungalCellRemoveIDs.push_back(id);
            setLastFungalCellChange(current_time);
            break;
        }
    }
}

void AgentManager::insertExistingComplexAtBoundary(Site *site, std::vector<Agent *> &complex, int initiatingCellId, double current_time) {

    Agent *phagocyte = nullptr;

    for (size_t i = 0; i < complex.size(); i++) {
        if (abm::util::isSubstring("ImmuneCell", complex[i]->getTypeName())) {
            phagocyte = complex[i];
        }
    }

    std::string agentType = phagocyte->getTypeName();

    const auto interactions = phagocyte->getInteractions()->getAllInteractions();
    int aliveCandida = phagocyte->getFeatureValueByName("aliveCandida");
    int killedCandida = phagocyte->getFeatureValueByName("killedCandida");
    int timestepOfFirstPhagocytosis = phagocyte->getFeatureValueByName("timestepOfFirstPhagocytosis");
    double radius = phagocyte->getMorphology()->getBasicSphereOfThis()->getRadius();
    //double interactionRadius = phagocyte->getInteractionRadius();


    Agent *positionTestAgent = 0;

    auto boundary_condition = site->get_boundary_condition();
    if (boundary_condition == "PeriodicBoundaries"){
        if (positionTestAgent != 0) {
            removeAgent(site, positionTestAgent, current_time);
        }
        Coordinate3D newPos = site->get_boundary()->boundary_cross_single_agent(phagocyte);
        //positionTestAgent = createAgent(site, agentType, newPos, current_time);
        positionTestAgent = emplace_back(site->getCellFactory()->createCell(agentType, std::make_unique<Coordinate3D>(newPos), generateNewID(), site, time_delta_, current_time, "")).get();
        if (positionTestAgent->getInteractions()->hasCollisions()){ // if complex has collisions in new position then do not move
            if (positionTestAgent != 0) {
                removeAgent(site, positionTestAgent, current_time);
            }
            positionTestAgent = emplace_back(site->getCellFactory()->createCell(agentType, std::make_unique<Coordinate3D>(phagocyte->getPreviousPosition()), generateNewID(), site, time_delta_, current_time, "")).get();
            //DEBUG_STDOUT("Shift of the complex not possible at the boundary");
        }
    }
    if (boundary_condition == "AperiodicBoundaries"){
        int rejections = 0;
        do {
            if (positionTestAgent != 0) {
                removeAgent(site, positionTestAgent, current_time);
            }
            Coordinate3D newPos = site->getRandomBoundaryPoint();
            if (newPos.x - radius <= site->getLowerLimits().x) {
                double correction = abs(newPos.x) + radius - abs(site->getLowerLimits().x);
                newPos = Coordinate3D{newPos.x + correction, newPos.y, newPos.z};
            }
            if (newPos.x + radius >= site->getUpperLimits().x) {
                double correction = abs(newPos.x) + radius - abs(site->getUpperLimits().x);
                newPos = Coordinate3D{newPos.x - correction, newPos.y, newPos.z};
            }
            if (newPos.y - radius <= site->getLowerLimits().y) {
                double correction = abs(newPos.y) + radius - abs(site->getLowerLimits().y);
                newPos = Coordinate3D{newPos.x, newPos.y + correction, newPos.z};
            }
            if (newPos.y + radius >= site->getUpperLimits().y) {
                double correction = abs(newPos.y) + radius - abs(site->getUpperLimits().y);
                newPos = Coordinate3D{newPos.x, newPos.y - correction, newPos.z};
            }
            if (newPos.z - radius <= site->getLowerLimits().z) {
                double correction = abs(newPos.z) + radius - abs(site->getLowerLimits().z);
                newPos = Coordinate3D{newPos.x, newPos.y, newPos.z + correction};
            }
            if (newPos.z + radius >= site->getUpperLimits().z) {
                double correction = abs(newPos.z) + radius - abs(site->getUpperLimits().z);
                newPos = Coordinate3D{newPos.x, newPos.y, newPos.z - correction};
            }
            Coordinate3D initialVector = site->getBoundaryInputVector();
            positionTestAgent = createAgent(site, agentType, newPos, &initialVector, current_time);


            rejections++;
        } while ((positionTestAgent->getInteractions()->hasCollisions() && rejections < 10000));
        if (rejections > 9999) {
            DEBUG_STDOUT("complex boundary cross: delete and new");
            positionTestAgent->getInteractions()->hasCollisions();
            DEBUG_STDOUT("too many agent-input-rejections on insert at boundary!");
            for (auto& i : complex) {
                DEBUG_STDOUT("at " << i->getPosition().x << " dim:" + std::to_string(site->getNumberOfSpatialDimensions()) + " site:" + site->getType());
                if (i != 0) {
                    removeAgent(site, i, current_time);
                }
            }
        }
    }
    else{
        DEBUG_STDOUT("NO BOUNDARY CONDITION TYPE INDICATED - can not handle boundary cross");
    }


    Coordinate3D newPos = positionTestAgent->getPosition();
    auto initCoord = std::make_unique<Coordinate3D>();
    *initCoord = newPos;
    auto newPhagocyte = site->getCellFactory()->createCell(agentType,
                                                           std::move(initCoord),
                                                           idHandling,
                                                           site,
                                                           time_delta_,
                                                           current_time, "");
    if (phagocyte->getId() != initiatingCellId) {
        site->getNeighbourhoodLocator()->removeSphereRepresentation(phagocyte->getMorphology()->getAllSpheresOfThis().front());
        site->getNeighbourhoodLocator()->updateDataStructures(newPhagocyte->getMorphology()->getAllSpheresOfThis().front());
    }
    //        cout << "[AgentManager] new phagocyte "<< newPhagocyte->getId() << '\n';
    if (newPhagocyte != nullptr) {
        allAgents.push_back(newPhagocyte);
        idHandling++;
    }
    //        System::getNeighbourhoodLocator()->removeSphereRepresentation(positionTestAgent->getAgentProperties()->getMorphology()->getAllSpheresOfThis().front())
    removeAgent(site, positionTestAgent, current_time);
    newPhagocyte->setFeatureValueByName("aliveCandida", aliveCandida);
    newPhagocyte->setFeatureValueByName("killedCandida", killedCandida);
    newPhagocyte->setFeatureValueByName("timestepOfFirstPhagocytosis", timestepOfFirstPhagocytosis);
    //        cout << "Phagocyte "<< agentType << " " << newPhagocyte->getId() << " killed Candida = " << newPhagocyte->getFeatureValueByName("killedCandida") << ", alive Candida = " << newPhagocyte->getFeatureValueByName("aliveCandida") << '\n';


    for (size_t i = 0; i < interactions.size(); i++) {
        if (interactions[i]->getInteractionName() == "PhagocyteFungusInteraction") {
            auto interaction = interactions[i];
            if (interaction->getCurrentState()->getStateName() == "Phagocytose"
                || interaction->getCurrentState()->getStateName() == "Lysis") {
                Cell *oldCandida = 0;
                if (phagocyte->getTypeName() == interaction->getFirstCell()->getTypeName()) {
                    oldCandida = interaction->getSecondCell();
                } else {
                    oldCandida = interaction->getFirstCell();
                }

                //old pathogen properties
                Coordinate3D nPos;
                std::string stateName = oldCandida->getCurrentCellState()->getStateName();
                std::string pathogenType = oldCandida->getTypeName();

                auto newInteractions = newPhagocyte->getInteractions()->getAllInteractions();
                std::vector<Agent *> newPhagocyteList;
                for (size_t j = 0; j < newInteractions.size(); j++) {
                    if (newInteractions[j]->getInteractionName().compare("PhagocyteFungusInteraction") == 0) {
                        if (newInteractions[j]->getCurrentState()->getStateName().compare("Phagocytose") == 0
                            || newInteractions[j]->getCurrentState()->getStateName().compare("Lysis") == 0) {
                            if (newPhagocyte->getTypeName().compare(newInteractions[j]->getFirstCell()->getTypeName())
                                == 0) {
                                newPhagocyteList.push_back(newInteractions[j]->getSecondCell());
                            } else {
                                newPhagocyteList.push_back(newInteractions[j]->getFirstCell());
                            }
                        }
                    }
                }
                nPos = newPhagocyte->getCoordinateWithinAgent(oldCandida);
                //
                Coordinate3D initialVector = site->getBoundaryInputVector();
                auto initCoord = std::make_unique<Coordinate3D>();
                *initCoord = nPos;

                //create new pathogen and insert it in list for agents
                auto newPathogen = std::dynamic_pointer_cast<FungalCell>(site->getCellFactory()->createCell(pathogenType,
                                                                                                            std::move(initCoord),
                                                                                                            idHandling,
                                                                                                            site,
                                                                                                            time_delta_,
                                                                                                            current_time, ""));
                if (newPathogen != 0) {
                    allAgents.push_back(newPathogen);
                    idHandling++;
                }
                newPathogen->changeState(stateName);
                //                    cout << "[AgentManager] newPathogen "<< pathogenType << " " << newPathogen->getId() << " " << newPathogen->getCurrentCellState()->getStateName() << '\n';
                //
                auto newInteraction = std::make_shared<PhagocyteFungusInteraction>(interaction->getIdentifier(),
                                                                                   newPhagocyte.get(),
                                                                                   newPathogen.get(),
                                                                                   true,
                                                                                   time_delta_,
                                                                                   current_time);
                std::string curStateName = interaction->getCurrentState()->getStateName();
                if (interaction->getCurrentState()->getStateName() == "Phagocytose") {
                    newInteraction->setState("Phagocytose");
                } else if (interaction->getCurrentState()->getStateName() == "Lysis") {
                    newInteraction->setState("Lysis");
                }
                //newPathogen->setPhagocytosed();
                newPhagocyte->getInteractions()->addInteraction(newInteraction);
                newPathogen->getInteractions()->addInteraction(newInteraction);
                if (oldCandida->getId() != initiatingCellId) {
                    site->getNeighbourhoodLocator()->removeSphereRepresentation(oldCandida->getMorphology()->getAllSpheresOfThis().front());
                    site->getNeighbourhoodLocator()->updateDataStructures(newPathogen->getMorphology()->getAllSpheresOfThis().front());
                }
            }
        }
    }

    //delete all old agents
    for (size_t i = 0; i < complex.size(); i++) {
        site->getAgentManager()->replaceAgent(site, complex[i], 0, current_time);
    }
    for (size_t i = 0; i < interactions.size(); i++) {
        interactions[i]->setDelted();
    }
}

void AgentManager::insertExistingAgentAtBoundary(Site *site, std::string agentType, std::string stateName, double current_time) {

    Coordinate3D initialPosition;
    Coordinate3D initialVector;

    Agent *agent = 0;
    int rejections = 0;
    do {
        if (agent != 0) {
            removeAgent(site, agent, current_time);
        }

        initialPosition = site->getRandomBoundaryPoint();
        initialVector = site->getBoundaryInputVector();
        agent = createAgent(site, agentType, initialPosition, &initialVector, current_time);
        agent->changeState(stateName);

        rejections++;
    } while ((agent->getInteractions()->hasCollisions() && rejections < 10000));

    if (rejections > 9999) {
        DEBUG_STDOUT("single agent boundary cross: new position");
        DEBUG_STDOUT("too many agent-input-rejections on insert at boundary!");
        DEBUG_STDOUT(
            "at " << "Position: (" << agent->getPosition().x << ", " << agent->getPosition().y << ", " << agent->getPosition().z << ") " << " dim:" +
                                                                                                                                                std::to_string(site->getNumberOfSpatialDimensions()) +
                                                                                                                                                " site:" + site->getType());
        if (agent != 0) {
            removeAgent(site, agent, current_time);
        }
    }
}
