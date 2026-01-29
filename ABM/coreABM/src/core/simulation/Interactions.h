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

#ifndef CORE_SIMULATION_INTERACTIONS_H
#define CORE_SIMULATION_INTERACTIONS_H

#include <memory>
#include <algorithm>
#include <list>
#include <map>
#include <vector>

#include "core/basic/Coordinate3D.h"
#include "core/simulation/neighbourhood/NeighbourhoodLocator.h"


class Cell;
class Interaction;
class InSituMeasurements;

class Interactions {
public:
// Class for wrapping all the interactions for one agent into on data container. In addition, this class provides utility functions for these interactions.
    Interactions(Cell *cell, NeighbourhoodLocator *nhLocator);
    virtual ~Interactions();

    void doWholeProcess(double timestep, double current_time, InSituMeasurements *measurments);
    void appendCollisionsToExistingInteractions(std::vector<std::shared_ptr<Collision>> &collisions);
    void removeCollisionsOfExistingInteractions(std::vector<std::shared_ptr<Collision>> &collisions);
    bool appendCollisionsToExistingInteractions(std::shared_ptr<Collision> collisions);
    void addNewInteractions(std::vector<std::shared_ptr<Collision>> &collisions, double time_delta, double current_time);
    void doAvoidanceInteractions(std::vector<std::shared_ptr<Collision>> &collisions, double time_delta,
                                 double current_time);
    void addInteraction(std::shared_ptr<Interaction> interaction);
    void executeAllInteractions(double timestep, double current_time);
    void removeClosedInteractions();
    void removeInteraction(Interaction *interaction);
    void removeAllInteractions();
    bool hasInteractions();
    bool hasPhagFungInteraction();
    bool hasCollisions();
    void avoidNewInteractions(double time_delta, double current_time);
    void displayInteractions();
    const std::vector<std::shared_ptr<Interaction>> &getAllInteractions();
    const std::map<Cell *, std::shared_ptr<Interaction>> &getAllInteractionPartners();

private:
    Cell *cell;
    NeighbourhoodLocator *neighbourhoodLocator;
    std::vector<std::shared_ptr<Interaction>> interactions;
    std::map<Cell *, std::shared_ptr<Interaction>> interactionPartners;
};

#endif /* CORE_SIMULATION_INTERACTIONS_H */
