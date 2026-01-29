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

#ifndef CORE_SIMULATION_ABSORBINGBOUNDARIES_H
#define CORE_SIMULATION_ABSORBINGBOUNDARIES_H

#include <memory>

#include "core/simulation/neighbourhood/BoundaryCondition.h"

class Analyser;

class AbsorbingBoundaries : public BoundaryCondition {
public:
  // Class for absorbing boundary conditions. This means that cells are removed from the system after they leave the environment. New cells are only inserted according to site-specific event.
    AbsorbingBoundaries(Site *site);

    void handleBoundaryCross(Agent *agent, Coordinate3D *moveVec, double current_time);
    std::string getTypeName();

    Coordinate3D update_pos(Coordinate3D pos) final;

};

#endif /* CORE_SIMULATION_ABSORBINGBOUNDARIES_H */
