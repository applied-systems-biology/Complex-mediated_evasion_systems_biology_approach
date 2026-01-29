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

#ifndef CORE_SIMULATION_BOUNDARYCONDITION_H
#define CORE_SIMULATION_BOUNDARYCONDITION_H

#include <memory>

#include "core/simulation/Agent.h"
#include "core/basic/Coordinate3D.h"

class Analyser;
class Site;

class BoundaryCondition {
public:
  // Abstract class for boundary condition
    explicit BoundaryCondition(Site *site);

    virtual ~BoundaryCondition() = default;
    virtual void handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time);
    virtual std::string getTypeName();
    virtual Coordinate3D boundary_cross_single_agent(Agent* agent);
    virtual Coordinate3D update_pos(Coordinate3D pos);
protected:
    Site *site;
};

#endif /* CORE_SIMULATION_BOUNDARYCONDITION_H */
