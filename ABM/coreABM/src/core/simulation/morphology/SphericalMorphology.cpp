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

#include <list>

#include "SphericalMorphology.h"
#include "SphereRepresentation.h"
#include "core/simulation/Site.h"

SphericalMorphology::SphericalMorphology() : MorphologyElement() {
}

SphericalMorphology::SphericalMorphology(Morphology *refMorphology, std::shared_ptr<Coordinate3D> position,
                                         double radius, std::string description, double creation_time) : MorphologyElement(refMorphology, description) {
    position_ = position;
    radius_ = radius;
    description_ = description;
    creation_time_ = creation_time;

    generateSphereRepresentation();
}

void SphericalMorphology::generateSphereRepresentation() {
    auto sR = std::make_unique<SphereRepresentation>(position_, radius_, this, description_, creation_time_);
    sphere_rep_.push_back(std::move(sR));
}

std::string SphericalMorphology::generatePovObject() {
    std::ostringstream ss;
    ss << "sphere {" << '\n';
    ss << "\t" << position_->x << "," << position_->y << "," << position_->z << "," << '\n';
    ss << "\t" << getRadius() << '\n';
    ss << "\t" << "pigment { color " << morphology_this_belongs_to_->getColorRGB()->printPovColorRGBT() << " }" << '\n';
    ss << "\t}" << '\n';

    return ss.str();
}