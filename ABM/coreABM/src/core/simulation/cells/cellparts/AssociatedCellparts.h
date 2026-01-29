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

#ifndef CORE_SIMULATION_ASSOCIATEDCELLPARTS_H
#define CORE_SIMULATION_ASSOCIATEDCELLPARTS_H

#include "core/simulation/Cell.h"
#include "core/simulation/morphology/Morphology.h"
#include "core/basic/Coordinate3D.h"

class AssociatedCellparts {
public:
    AssociatedCellparts() =default;
    
    AssociatedCellparts(Cell *mothercell, std::shared_ptr<Coordinate3D> connectionPoint);

    virtual ~AssociatedCellparts() = default;
   
    Coordinate3D* getMcConnectionPoint();
    
protected:

    Cell *mothercell;
    std::shared_ptr<Coordinate3D> mcConnectionPoint;
    
};

#endif /* CORE_SIMULATION_ASSOCIATEDCELLPARTS_H */
