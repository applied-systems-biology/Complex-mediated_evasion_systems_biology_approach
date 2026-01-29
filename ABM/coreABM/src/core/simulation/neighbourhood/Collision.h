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

#ifndef CORE_SIMULATION_COLLISION_H
#define CORE_SIMULATION_COLLISION_H

#include <map>
#include <set>

#include "core/simulation/morphology/SphereRepresentation.h"

class Cell;

enum class MeasurementType {
    NO_INTERACTION,
    EXISTING_INTERACTION,
    NEW_INTERACTION
};

class Collision {
public:
  // Class for collision check between cells represented as spheres
    Collision() = default;
    Collision(Cell *collisionCell, SphereRepresentation *mySphere, SphereRepresentation *collisionSphere,
              double overlap);
    [[nodiscard]] Cell *getCell() const { return cell; };
    [[nodiscard]] Cell *getCollisionCell() const { return collisionCell; };
    SphereRepresentation *getMySphere() { return mySphere; };
    SphereRepresentation *getCollisionSphere() { return collisionSphere; };
    double getTimeToFirstContact() { return timeToFirstContact; };
    void calculateTimeTillFirstContact();

    MeasurementType type{MeasurementType::NO_INTERACTION};
private:
    Cell *cell{};
    Cell *collisionCell{};
    SphereRepresentation *mySphere{};
    SphereRepresentation *collisionSphere{};
    double timeToFirstContact{};
    double overlap{};
};

#endif /* CORE_SIMULATION_COLLISION_H */
