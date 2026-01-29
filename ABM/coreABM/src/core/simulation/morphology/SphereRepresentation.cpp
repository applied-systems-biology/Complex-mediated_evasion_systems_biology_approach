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

#include "SphereRepresentation.h"
#include "MorphologyElement.h"
#include "core/simulation/Cell.h"


SphereRepresentation::SphereRepresentation() {
}

SphereRepresentation::SphereRepresentation(std::shared_ptr<Coordinate3D> coord, double radius,
                                           MorphologyElement *morphologyElement, std::string description, double creation_time) {
    this->position = coord;
    this->morphologyElementThisBelongsTo = morphologyElement;
    this->radius = radius;
    this->description_ = description;
    this->creation_time_ = creation_time;

    id = morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->getAgentManager()->getNextSphereRepresentationId(this);

    morphologyElementThisBelongsTo->getMorphologyThisBelongsTo()->getCellThisBelongsTo()->getSite()->getNeighbourhoodLocator()->addSphereRepresentation(this);

}

SphereRepresentation::~SphereRepresentation() {
}

Coordinate3D SphereRepresentation::getEffectiveConnection(SphereRepresentation *sphereRep) {
    return Coordinate3D(sphereRep->getPosition() - *position);
}

void SphereRepresentation::shiftPosition(Coordinate3D *shifter) {
    *position += *shifter;
}

void SphereRepresentation::setRadiusToOrigin(double r) {
    SphericCoordinate3D oldPos = abm::util::toSphericCoordinates(*position);
    oldPos.r = r;
    *position = abm::util::toCartesianCoordinates(oldPos);
}
