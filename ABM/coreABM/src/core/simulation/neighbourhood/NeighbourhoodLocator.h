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

#ifndef CORE_SIMULATION_NEIGHBOURHOODLOCATOR_H
#define CORE_SIMULATION_NEIGHBOURHOODLOCATOR_H

#include <set>
#include <algorithm>
#include <string>
#include <vector>

#include "core/simulation/morphology/SphereRepresentation.h"


class Agent;
class Collision;
class Site;

class NeighbourhoodLocator {
public:
    /// Abstract class for a cell neighbourhood locator to speed up collision checks between cells
    NeighbourhoodLocator(Site *Site);
    virtual ~NeighbourhoodLocator();
    virtual void instantiate();
    virtual std::vector<std::shared_ptr<Collision>> getCollisions(Agent *agent);
    virtual bool hasCollision(Agent *agent) { return false; };
    virtual std::vector<Coordinate3D> getCollisionSpheres(SphereRepresentation *sphereRep, Coordinate3D dirVec);
    virtual void updateDataStructures(SphereRepresentation *sphereRep);
    virtual void removeSphereRepresentation(SphereRepresentation *sphereRep);
    virtual void addSphereRepresentation(SphereRepresentation *sphereRep);
    virtual int controlFunction() { return 0; };
    virtual std::string getTypeName();
    virtual void check() {};
    virtual int getNumberOfAgentTypeInBalloonList(std::string agentType) { return 0; };
protected:
    Site *site_;
};

#endif /* CORE_SIMULATION_NEIGHBOURHOODLOCATOR_H */
