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

#ifndef CORE_SIMULATION_AGENT_H
#define CORE_SIMULATION_AGENT_H

#include <iostream>

#include "core/basic/Coordinate3D.h"
#include "core/simulation/movement/RandomWalk.h"
#include "core/simulation/movement/Movement.h"
#include "core/simulation/morphology/Morphology.h"
#include "core/simulation/states/CellState.h"

class Site; //forward declaration
class Interactions; //forward declaration

class Agent {
public:
  // Abstract class for agents in the hABM. This class provides the cell class with it main functionality.
    Agent();
    Agent(std::unique_ptr<Coordinate3D>, int, Site *);
    virtual ~Agent() = default;

    void setPosition(Coordinate3D newPos);
    void setBeenMovedThisTimestep(bool newHasBeenMoved);
    void setTimestepLastTreatment(double current_time);
    void setDeleted();
    void setPreviousPosition(Coordinate3D *pos) { *previousPosition = *pos; }
    void setInputRate(double input_rate) {inputRate = input_rate;}
    void setMorphology(std::shared_ptr<Morphology> morphology) {morphology_ = std::move(morphology);};
    void setMovement(std::shared_ptr<Movement> movement) {movement_ = std::move(movement);};
    void setPassiveMovement(std::shared_ptr<Movement> passiveMovement) {passive_movement_ = std::move(passiveMovement);};
    void setInteractions(std::shared_ptr<Interactions> interactions) {interactions_ = interactions;};
    void resetAgent(Coordinate3D, double current_time);

    [[nodiscard]] int getId() const;
    [[nodiscard]] bool isDeleted() const { return is_deleted_; }
    [[nodiscard]] bool agentTreatedInCurrentTimestep(double current_time) const;
    bool hasBeenMovedThisTimestep();
    bool coordinateIsInsideAgent(Coordinate3D *, Agent *);
    bool hasCollisionsInsideAgent(Agent *, const std::vector<Agent *> &);
    bool shiftPosition(Coordinate3D *shifter,
                       double current_time,
                       SphereRepresentation *sphereRep = 0,
                       std::string originCall = "not known");
    double getLifetime(double);
    double getInitialTime();
    double getInputRate() {return inputRate;};
    Site *getSite();
    Coordinate3D getPosition();
    Coordinate3D *getCurrentShift();
    Coordinate3D getPreviousPosition() { return *previousPosition; };
    Coordinate3D getInitialPosition() { return *initialPosition; };
    Coordinate3D getCurrentPosition() { return *position; };
    Coordinate3D getCoordinateWithinAgent(Agent *);
    Movement *getMovement() {return movement_.get();};
    Movement *getPassiveMovement() {return passive_movement_.get();};
    Interactions *getInteractions() {return interactions_.get();};
    Morphology *getMorphology() {return morphology_.get();};

    std::map<std::string, double> molecule_uptake;

    virtual void doAllActionsForTimestep(double timestep, double current_time) = 0;
    virtual void move(double timestep, double current_time) = 0;
    virtual std::string getTypeName() = 0;
    virtual CellState *getCurrentCellState() = 0;
    virtual std::string generatePovObject() = 0;
    virtual void setPassive() = 0;
    virtual void setVariableOnEvent(std::string variable, double value) = 0;
    virtual void setFeatureValueByName(std::string featureName, int value) = 0;
    virtual void applyMethodByName(std::string mehtodName) = 0;
    virtual void changeState(std::string stateName) = 0;
    virtual double getFeatureValueByName(std::string featureName) = 0;
    virtual std::shared_ptr<CellState> getCellStateByName(std::string nameOfState) = 0;
    virtual void setState(std::shared_ptr<CellState> state) = 0;
    virtual Morphology *getSurface() = 0;
    virtual Coordinate3D get_gradient() = 0;
    std::vector<Agent*> getComplex();

protected:
    void setInitialPosition(Coordinate3D initPos);
    unsigned int id{};
    bool positionShiftAllowed{};
    bool hasBeenMoved{};
    bool passive{};
    bool is_deleted_{};
    double initialTime{};
    double timestepLastTreatment{};
    double inputRate{};

    std::shared_ptr<Movement> movement_{};
    std::shared_ptr<Movement> passive_movement_{};
    std::shared_ptr<Interactions> interactions_{};
    std::shared_ptr<Morphology> morphology_{};

    Site *site;
    std::unique_ptr<Coordinate3D> currShift{};
    std::unique_ptr<Coordinate3D> initialPosition{};
    std::unique_ptr<Coordinate3D> previousPosition{};
    std::shared_ptr<Coordinate3D> position{};
    std::shared_ptr<Movement> movement{};
    std::shared_ptr<Movement> passiveMovement{};
};

#endif /* CORE_SIMULATION_AGENT_H */
