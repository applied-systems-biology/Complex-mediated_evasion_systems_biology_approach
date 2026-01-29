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

#ifndef CORE_SIMULATION_SPHEREREPRESENTATION_H
#define CORE_SIMULATION_SPHEREREPRESENTATION_H

#include <memory>
#include <sstream>

#include "core/basic/Coordinate3D.h"
#include "core/basic/ColorRGB.h"

class MorphologyElement;

class SphereRepresentation {
public:
  // Class for a spherical representation based on a spherical morphology.
    SphereRepresentation();
    SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius, MorphologyElement *morphologyElement,
                         std::string description = "", double creation_time=0.0);

    virtual ~SphereRepresentation();
    Coordinate3D getPosition() { return *position; };
    double getRadius() { return radius; };
    double getCreationTime() { return creation_time_; }
    void setRadius(double r) { radius = r; };
    int getId() { return id; };
    std::string getDescription() { return description_; };
    MorphologyElement *getMorphologyElementThisBelongsTo() { return morphologyElementThisBelongsTo; };
    Coordinate3D getEffectiveConnection(SphereRepresentation *sphereRep);
    void shiftPosition(Coordinate3D *shifter);
    void setRadiusToOrigin(double r);

private:
    std::shared_ptr<Coordinate3D> position;
    double radius;
    double creation_time_;
    std::string description_;
    MorphologyElement *morphologyElementThisBelongsTo;
    int id;
};

#endif /* CORE_SIMULATION_SPHEREREPRESENTATION_H */
