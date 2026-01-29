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

/*
 * File:   PeriodicBoundaries.cpp
 * Author: jpollmae
 * 
 * Created on 27. Februar 2012, 14:51
 */

#include "PeriodicBoundaries.h"
#include "core/simulation/Site.h"
#include "core/simulation/AgentManager.h"
#include "core/simulation/Interactions.h"


PeriodicBoundaries::PeriodicBoundaries(Site* site) : BoundaryCondition(site) {
}



void PeriodicBoundaries::handleBoundaryCross(Agent* agent, Coordinate3D* moveVec, double current_time) {
    std::string agentType = agent->getTypeName();
    Agent* test_agent = 0;
    auto new_coord = boundary_cross_single_agent(agent);
    auto newCoord = std::make_unique<Coordinate3D>(new_coord);
    site->getAgentManager()->replaceAgent(site, agent, std::move(newCoord), moveVec, current_time);
}
Coordinate3D PeriodicBoundaries::boundary_cross_single_agent(Agent* agent){
    Coordinate3D agentsCoord = agent->getPosition();
    // cout << "old:" << agentsCoord->printCoordinate() << '\n';
    double newX, newY, newZ, diffX, diffY, diffZ;
    bool lower, upper;

    Coordinate3D pointA = site->getLowerLimits();
    Coordinate3D pointB = site->getUpperLimits();
      switch (site->getNumberOfSpatialDimensions()) {
      case 2:
          if (!((lower = (agentsCoord.x >= pointA.x)) && (upper = (agentsCoord.x <= pointB.x)))) {
              if (!lower) {
                  diffX = agentsCoord.x - pointA.x;
                  newX = pointB.x + diffX;
              } else {
                  diffX = agentsCoord.x - pointB.x;
                  newX = pointA.x + diffX;
              }
          } else {
              newX = agentsCoord.x;
          }

          if (!((lower = (agentsCoord.y >= pointA.y)) && (upper = (agentsCoord.y <= pointB.y)))) {
              if (!lower) {
                  diffY = agentsCoord.y - pointA.y;
                  newY = pointB.y + diffY;
              } else {
                  diffY = agentsCoord.y - pointB.y;
                  newY = pointA.y + diffY;
              }
          } else {
              newY = agentsCoord.y;
          }

          newZ = agentsCoord.z;
          break;

      case 3:
          if (!((lower = (agentsCoord.x >= pointA.x)) && (upper = (agentsCoord.x <= pointB.x)))) {
              if (!lower) {
                  diffX = agentsCoord.x - pointA.x;
                  newX = pointB.x + diffX;
              } else {
                  diffX = agentsCoord.x - pointB.x;
                  newX = pointA.x + diffX;
              }
          } else {
              newX = agentsCoord.x;
          }

          if (!((lower = (agentsCoord.y >= pointA.y)) && (upper = (agentsCoord.y <= pointB.y)))) {
              if (!lower) {
                  diffY = agentsCoord.y - pointA.y;
                  newY = pointB.y + diffY;
              } else {
                  diffY = agentsCoord.y - pointB.y;
                  newY = pointA.y + diffY;
              }
          } else {
              newY = agentsCoord.y;
          }

          if (!((lower = (agentsCoord.z >= pointA.z)) && (upper = (agentsCoord.z <= pointB.z)))) {
              if (!lower) {
                  diffZ = agentsCoord.z - pointA.z;
                  newZ = pointB.z + diffZ;
              } else {
                  diffZ = agentsCoord.z - pointB.z;
                  newZ = pointA.z + diffZ;
              }
          } else {
              newZ = agentsCoord.z;
          }
          double radius;
          radius = agent->getMorphology()->getBasicSphereOfThis()->getRadius();
          if (newX - radius <= site->getLowerLimits().x) {
              double correction = abs(newX) + radius - abs(site->getLowerLimits().x);
              newX += correction;
          }
          if (newX + radius >= site->getUpperLimits().x) {
              double correction = abs(newX) + radius - abs(site->getUpperLimits().x);
              newX -= correction;
          }
          if (newY - radius <= site->getLowerLimits().y) {
              double correction = abs(newY) + radius - abs(site->getLowerLimits().y);
              newY += correction;
          }
          if (newY + radius >= site->getUpperLimits().y) {
              double correction = abs(newY) + radius - abs(site->getUpperLimits().y);
              newY -= correction;
          }
          if (newZ - radius <= site->getLowerLimits().z) {
              double correction = abs(newZ) + radius - abs(site->getLowerLimits().z);
              newZ += correction;
          }
          if (newZ + radius >= site->getUpperLimits().z) {
              double correction = abs(newZ) + radius - abs(site->getUpperLimits().z);
              newZ -= correction;
          }
          break;

      default:
          if (!((lower = (agentsCoord.x >= pointA.x)) && (upper = (agentsCoord.x <= pointB.x)))) {
              if (!lower) {
                  diffX = agentsCoord.x - pointA.x;
                  newX = pointB.x + diffX;
              } else {
                  diffX = agentsCoord.x - pointB.x;
                  newX = pointA.x + diffX;
              }
          } else {
              newX = agentsCoord.x;
          }

          if (!((lower = (agentsCoord.y >= pointA.y)) && (upper = (agentsCoord.y <= pointB.y)))) {
              if (!lower) {
                  diffY = agentsCoord.y - pointA.y;
                  newY = pointB.y + diffY;
              } else {
                  diffY = agentsCoord.y - pointB.y;
                  newY = pointA.y + diffY;
              }
          } else {
              newY = agentsCoord.y;
          }

          if (!((lower = (agentsCoord.z >= pointA.z)) && (upper = (agentsCoord.z <= pointB.z)))) {
              if (!lower) {
                  diffZ = agentsCoord.z - pointA.z;
                  newZ = pointB.z + diffZ;
              } else {
                  diffZ = agentsCoord.z - pointB.z;
                  newZ = pointA.z + diffZ;
              }
          } else {
              newZ = agentsCoord.z;
          }
          break;
      }

      // cout << agent->getPosition()->printCoordinate() << "" << agent->getCurrentShift()->printVector() << " " << moveVec->printVector() << '\n';
      return Coordinate3D{newX, newY, newZ};
}

std::string PeriodicBoundaries::getTypeName(){
  return "PeriodicBoundaries";
}
Coordinate3D PeriodicBoundaries::update_pos(Coordinate3D pos) {
    auto upper_limit = site->getUpperLimits();
    auto lower_limit = site->getLowerLimits();

    if (pos.x < lower_limit.x){
        pos.x += upper_limit.x - lower_limit.x;
    }
    else if (pos.x > upper_limit.x){
        pos.x -= upper_limit.x - lower_limit.x;
    }

    if (pos.y < lower_limit.y){
        pos.y += upper_limit.y - lower_limit.y;
    }
    else if (pos.y > upper_limit.y){
        pos.y -= upper_limit.y - lower_limit.y;
    }

    if (pos.z < lower_limit.z){
        pos.z += upper_limit.z - lower_limit.z;
    }
    else if (pos.z > upper_limit.z){
        pos.z -= upper_limit.z - lower_limit.z;
    }
    return pos;
}
