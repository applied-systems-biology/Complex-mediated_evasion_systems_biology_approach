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

#ifndef APERIODICBOUNDARIES_H
#define	APERIODICBOUNDARIES_H

#include "core/simulation/neighbourhood/BoundaryCondition.h"

class AperiodicBoundaries : public BoundaryCondition{
public:
    AperiodicBoundaries(Site* site);
    void handleBoundaryCross(Agent* agent, Coordinate3D* moveVec, double current_time);

  std::string getTypeName();

  Coordinate3D update_pos(Coordinate3D pos) final;

private:

    /**
     * transforms a given previous movement according to a randomly chosen boundary position
     * that the speed but not the direction of the speed is conserved
     * @param prevMove previous movement of this agent 
     * @param boundaryPoint a given randomly chosen boundary point
     * @return a new movement with conserved speed of the previous move
     */
    Coordinate3D aperiodicMovementTransformation(Coordinate3D* prevMove, Coordinate3D* boundaryPoint, double delta_time);
    
};

#endif	/* APERIODICBOUNDARIES_H */

