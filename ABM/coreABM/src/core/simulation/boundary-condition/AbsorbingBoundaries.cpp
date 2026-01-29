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

#include "AbsorbingBoundaries.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Site.h"
#include "core/analyser/InSituMeasurements.h"

AbsorbingBoundaries::AbsorbingBoundaries(Site *site)
        : BoundaryCondition(site) {}

void AbsorbingBoundaries::handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time) {
    site->getAgentManager()->replaceAgent(site, agent, nullptr, current_time);
}

std::string AbsorbingBoundaries::getTypeName() {
    return "AbsorbingBoundaries";
}

Coordinate3D AbsorbingBoundaries::update_pos(Coordinate3D pos) {
}
