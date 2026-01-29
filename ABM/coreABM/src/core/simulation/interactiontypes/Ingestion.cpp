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

#include "Ingestion.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Site.h"
#include "core/simulation/Interactions.h"


void Ingestion::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {
    //active cell is the ingested one (the one that is eaten)
    //passive cell is the eater-cell

    Cell *passiveCell, *activeCell;

    if (isIngestingType(cell)) {
        activeCell = interaction->getOtherCell(cell);
        passiveCell = cell;
    } else {
        activeCell = cell;
        passiveCell = interaction->getOtherCell(cell);
    }

    Coordinate3D shift = passiveCell->getEffectiveConnection(activeCell);
    double currentDistance = shift.getMagnitude();

    // put the ingested cell at r/2 distance of the center of the ingesting cell
    double multiplicatorForNormalVec =
            passiveCell->getSurface()->getAllSpheresOfThis().front()->getRadius() / (2 * currentDistance) - 1;
    shift *= multiplicatorForNormalVec;
    activeCell->shiftPosition(&shift, current_time, 0, getTypeName());

    // prevent that ingested cell goes out of the site
    Site *currentSite = cell->getSite();
    std::string siteName = currentSite->getType();

    // Ensure that passive cell does not get out of spheric site
//    auto upperCheck = currentSite->getUpperLimits() - activeCell->getPosition();
//    auto lowerCheck = activeCell->getPosition() - currentSite->getLowerLimits();
//
//


    std::shared_ptr<Collision> currentCollision;
    while ((currentCollision = interaction->getNextCollision()) != 0) {
        currentCollision = nullptr;
    }
}

bool Ingestion::isIngestingType(Cell *cell) {
    bool result = false;

    if (abm::util::isSubstring("ImmuneCell", cell->getTypeName())) { result = true; }

    return result;
}

std::string Ingestion::getTypeName() const {
    return "Ingestion";
}