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

#include "RigidContacting.h"
#include "core/simulation/Interaction.h"


void RigidContacting::handleInteraction(Interaction *interaction, Cell *cell, double timestep, double current_time) {

    std::shared_ptr<Collision> currentCollision = nullptr;
    std::shared_ptr<Collision> collisionToHandle = nullptr;

    //check for the collision with the lowest time til first collision
    while ((currentCollision = interaction->getNextCollision()) != 0) {
        if (collisionToHandle != 0) {
            if (currentCollision->getTimeToFirstContact() < collisionToHandle->getTimeToFirstContact()) {
                collisionToHandle = currentCollision;
            } else {
                currentCollision = nullptr;
            }
        } else {
            collisionToHandle = currentCollision;
        }
    }
    double timeStepFirstContact = collisionToHandle->getTimeToFirstContact();

    //go back the current movement path until only the first contact between the spheres occurs
    Coordinate3D back_shift = *cell->getMovement()->getCurrentMove();
    back_shift *= -1.0 * (timestep - timeStepFirstContact) / timestep;

    //desiredShift->addVector(&backShift);
    cell->shiftPosition(&back_shift, current_time, 0, getTypeName());

    collisionToHandle = nullptr;
}