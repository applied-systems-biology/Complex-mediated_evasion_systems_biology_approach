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

#ifndef CORE_SIMULATION_CELL_H
#define CORE_SIMULATION_CELL_H

#include <memory>

#include "core/simulation/Agent.h"
#include "core/utils/io_util.h"
#include "core/simulation/Site.h"

class Analyser;

class Cell : public Agent {
public:
    /// Class for an abstract cell that is further specified by inherited classes
    Cell(std::unique_ptr<Coordinate3D>, int, Site *, double time_delta, double current_time);

    /*!
     * Central function: Performs all actions for one timestep for one cell
     * @param timestep Double for courrent timestep
     * @param current_time Double for current time
     */
    void doAllActionsForTimestep(double timestep, double current_time);

    void setState(std::shared_ptr<CellState> cstate) final;
    void setPassive() override { passive = true; }
    void move(double timestep, double current_time) override;
    void setFeatureValueByName(std::string featureName, int value) override;
    void setVariableOnEvent(std::string variable, double value);
    void applyMethodByName(std::string mehtodName) override;
    void changeState(std::string stateName) override;
    void recieveInteractionEvent(InteractionEvent *ievent, double current_time);
    void setInitialState(double time_delta, double current_time);
    void setExistingState(std::string stateName, double time_delta, double current_time);
    void addIngestions(int id);

    double getFeatureValueByName(std::string featureName);
    Coordinate3D get_gradient() final { return cumulative_persistence_gradient; }
    Coordinate3D getEffectiveConnection(Cell *cell);
    std::string getTypeName() override;
    std::string generatePovObject() final;
    CellState *getCurrentCellState() final;
    std::shared_ptr<CellState> getCellStateByName(std::string nameOfState) final;
    Morphology *getSurface();
    Interactions *getInteractions();
    virtual int get_nb_bound(){return 0;}

    virtual void handleControlledAgents(double timestep);
    virtual void passiveMove(double timestep, double current_time);
    virtual void doMorphologicalChanges(double timestep, double current_time){};
    virtual void interactWithMolecules(double timestep){};
    virtual void setup(double time_delta,
                       double current_time,
                       abm::util::SimulationParameters::AgentParameters *parameters);
    virtual abm::util::SimulationParameters::FungalParameters* getFungalParameter() {return nullptr;};


protected:
    virtual void handleInteractionEvent(InteractionEvent *ievent, double current_time);
    std::vector<int> ingestionCounter{};
    std::shared_ptr<Morphology> surface{};
    std::shared_ptr<Interactions> interactions{};
    std::map<std::string, std::shared_ptr<CellState>> cellStates{};
    std::shared_ptr<CellState> cellState{};
    Coordinate3D cumulative_persistence_gradient{};

};

#endif    /* CORE_SIMULATION_CELL_H */

