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

#include "BoundaryCondition.h"
#include "core/simulation/Site.h"
#include "core/analyser/Analyser.h"


BoundaryCondition::BoundaryCondition(Site *site) {
    this->site = site;
}

void BoundaryCondition::handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time) {
}

std::string BoundaryCondition::getTypeName() {
    return "BoundaryCondition";
}

Coordinate3D BoundaryCondition::boundary_cross_single_agent(Agent* agent){
}
Coordinate3D BoundaryCondition::update_pos(Coordinate3D pos) {
}
