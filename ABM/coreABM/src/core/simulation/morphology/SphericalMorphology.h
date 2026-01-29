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

#ifndef CORE_SIMULATION_SPHERICALMORPHOLOGY_H
#define CORE_SIMULATION_SPHERICALMORPHOLOGY_H

#include "MorphologyElement.h"
#include "SphereRepresentation.h"

class SphericalMorphology : public MorphologyElement {
public:
  // Class for a spherical morphology
    SphericalMorphology();
    SphericalMorphology(Morphology *refMorphology, std::shared_ptr<Coordinate3D> position, double radius,
                        std::string description = "", double creation_time=0.0);

    void generateSphereRepresentation() final;

    std::string generatePovObject();
    double getRadius() { return sphere_rep_.front()->getRadius(); };

    std::string getTypeName() { return "SphericalMorphology"; };

private:
    double radius_;
};

#endif /* CORE_SIMULATION_SPHERICALMORPHOLOGY_H */
