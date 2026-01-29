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

#ifndef CORE_SIMULATION_MORPHOLOGY_H
#define CORE_SIMULATION_MORPHOLOGY_H

#include <list>

#include "core/basic/ColorRGB.h"
#include "core/basic/Coordinate3D.h"
#include "core/simulation/morphology/MorphologyElement.h"

class Cell;

class Morphology {
public:
  // Class for handling the morphology of cells
    Morphology() = default;

    Morphology(const std::string &color, Cell *cell);
    virtual ~Morphology() = default;
    virtual std::string getTypeName();
    std::string generatePovObject();
    double getMaximumDistanceFromCenter();
    ColorRGB *getColorRGB();
    void setColorRGB(std::unique_ptr<ColorRGB> col);
    Cell *getCellThisBelongsTo() { return cell_this_belongs_to_; };
    void appendAssociatedCellpart(std::unique_ptr<MorphologyElement> morphElement);
    std::vector<SphereRepresentation *> getAllSpheresOfThis();
    SphereRepresentation *getBasicSphereOfThis();
    double getVolume();

protected:
    std::unique_ptr<ColorRGB> color_rgb_{};
    std::list<std::unique_ptr<MorphologyElement>> morphologyElements;
    Cell *cell_this_belongs_to_{};

};

#endif /* CORE_SIMULATION_MORPHOLOGY_H */
