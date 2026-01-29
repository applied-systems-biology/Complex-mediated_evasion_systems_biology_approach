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

#ifndef CORE_SIMULATION_MORPHOLOGYELEMENT_H
#define CORE_SIMULATION_MORPHOLOGYELEMENT_H

#include <utility>
#include <vector>

#include "core/basic/Coordinate3D.h"
#include "core/simulation/morphology/SphereRepresentation.h"

class Morphology;

class MorphologyElement {
public:
  // Abstract class for morphology representation
    MorphologyElement() = default;
    MorphologyElement(Morphology *morphologyRef, std::string description)
            : morphology_this_belongs_to_(morphologyRef), description_(std::move(description)) {};
    virtual ~MorphologyElement() = default;
    virtual void generateSphereRepresentation() = 0;

    const std::vector<std::shared_ptr<SphereRepresentation>> &getSphereRepresentation();
    std::string getDescription();
    Morphology *getMorphologyThisBelongsTo();
    double getVolume() const;

protected:
    std::shared_ptr<Coordinate3D> position_{};
    double creation_time_{};
    ColorRGB color;
    double volume_{};
    Morphology *morphology_this_belongs_to_{};
    std::string description_;
    std::vector<std::shared_ptr<SphereRepresentation>> sphere_rep_;
};

#endif /* CORE_SIMULATION_MORPHOLOGYELEMENT_H */
