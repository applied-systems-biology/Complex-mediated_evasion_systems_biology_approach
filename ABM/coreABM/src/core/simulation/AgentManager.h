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

#ifndef CORE_SIMULATION_AGENTMANAGER_H
#define CORE_SIMULATION_AGENTMANAGER_H

#include <map>
#include <memory>
#include <set>
#include <utility>
#include <vector>

#include "core/simulation/morphology/SphereRepresentation.h"
#include "core/utils/io_util.h"

class Site;
class Analyser;
class Agent;
class Cell;

class AgentManager {
public:
    /// Class for managing all the agents in data container. It provides functionality for input and replacement of agents during a simulation.
    AgentManager(abm::util::SimulationParameters& parameters, Site *site);

    using iterator = std::vector<std::shared_ptr<Agent>>::iterator;
    using const_iterator = std::vector<std::shared_ptr<Agent>>::const_iterator;
    iterator begin() { return allAgents.begin(); }
    iterator end() { return allAgents.end(); }
    [[nodiscard]] const_iterator begin() const { return allAgents.begin(); }
    [[nodiscard]] const_iterator end() const { return allAgents.end(); }
    std::shared_ptr<Agent> &emplace_back(std::shared_ptr<Agent> &&value) {return allAgents.emplace_back(std::forward<std::shared_ptr<Agent>>(value));}

    /*!
     * Exponentially distributed inputs of agents for a rate lamba which can be defined for each agent type
     * @param current_time Double that contains current time
     * @param random_generator Randomizer object that contains grandom generator
     */
    void inputOfAgents(double current_time, Randomizer *random_generator);

    /*!
     * Inserts an agent at the systems boundary (e.g. alveolar entrance ring or pores of kohn)
     * @param site Site object (e.g. AlveoleSite)
     * @param agentType String that contains agent type
     */
    void insertAgentAtBoundary(Site *site, std::string agentType, double current_time);

    /*!
     * Replaces an agent with a new agent at a given position and move
     * @param site Site object (e.g. AlveoleSite)
     * @param agentToReplace Agent object that contains agent that is replaced
     * @param c Coordinate3D object that contains coordinates of new agent
     * @param prevMove Coordinate3D object that contains corrdinates of previous move
     */
    void replaceAgent(Site *site, Agent *agentToReplace, std::unique_ptr<Coordinate3D> c, Coordinate3D * prevMove, double current_time);

    /*!
     * Replaces an agent with another agent
     * @param site Site object (e.g. AlveoleSite)
     * @param agentToReplace Agent object that contains agent that is replaced
     * @param newAgent Agent object that replaces the previous agent
     */
    void replaceAgent(Site *site, Agent *agentToReplace, std::shared_ptr<Agent> newAgent, double current_time);

    /// Removes agent from the system
    void removeAgent(Site *site, Agent *agent, double current_time);

    /// Removes agents from the system that were previously set to deleted
    void cleanUpAgents(double current_time);

    /// Removes a sphere object
    void removeSphereRepresentation(SphereRepresentation *sphereRep);
    void incrementIdHandling();
    void addFungalCellToList(Agent *fungus) { activeFungalCells.emplace_back(fungus); };
    virtual void removeFungalCellFromList(int id, double current_time);

    Agent *createAgent(Site *, std::string, Coordinate3D, Coordinate3D *, double current_time);
    int generateNewID() { return idHandling++;}

    void setLambdaInput(double speed, double persistenceTime);
    void setInitFungalQuantity() { initFungalQuantity = activeFungalCells.size(); }
    void setLastFungalCellChange(double lcc) { lastFungalCellChange = lcc; }
    int getAgentQuantity(std::string agenttype);
    int getNextSphereRepresentationId(SphereRepresentation *sphereRep);
    [[nodiscard]] double getLastFungalCellChange() const { return lastFungalCellChange; };
    [[nodiscard]] int getIdHandling() const;
    [[nodiscard]] int getInitFungalQuantity() const { return initFungalQuantity; }
    double getOccupancyDensityOfSpace();
    Cell *getCellBySphereRepId(int sphereRepId);
    const std::vector<std::shared_ptr<Agent>> &getAllAgents();
    SphereRepresentation *getSphereRepBySphereRepId(int sphereRepId);
    std::set<SphereRepresentation *> *getAllSphereRepresentations() { return &allSphereRepresentations; };
    std::vector<Agent *> getAllFungalCells() { return activeFungalCells; };
    std::vector<std::string> getAllAgentTypes();
    void insertExistingAgentAtBoundary(Site *site, std::string agentType, std::string stateName, double current_time);
    void insertExistingComplexAtBoundary(Site *site, std::vector<Agent *> &complex, int initiatingCellId, double current_time);

protected:
    std::vector<std::shared_ptr<Agent>> allAgents;
    std::map<int, Cell *> sphereIdToCell;
    std::map<int, SphereRepresentation *> sphereIdToSphereRep;
    std::set<SphereRepresentation *> allSphereRepresentations;
    int idHandling;
    int idHandlingSphereRepresentation;
    double lastFungalCellChange{};
    double lastInputEventTime{};
    double nextInputEventTime{};
    unsigned int lastQuantity{};
    std::vector<Agent *> activeFungalCells{};
    int initFungalQuantity{};
    double lambdaInput{};
    std::vector<double> fungalCellRemoveTimes{};
    std::vector<double> fungalCellRemoveIDs{};
    Site *site{};
    double time_delta_{};
    std::vector<std::string> agent_types_{};
};

#endif /* CORE_SIMULATION_AGENTMANAGER_H */
