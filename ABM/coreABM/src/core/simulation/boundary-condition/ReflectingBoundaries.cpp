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
 * File:   ReflectingBoundaries.cpp
 * Author: jpollmae
 * 
 * Created on 27. Februar 2012, 14:27
 */

#include "ReflectingBoundaries.h"

#include "core/basic/Coordinate3D.h"
#include "core/simulation/Site.h"

ReflectingBoundaries::ReflectingBoundaries(Site* site) : BoundaryCondition(site) {
}

void ReflectingBoundaries::handleBoundaryCross(Agent* agent, Coordinate3D* moveVec, double current_time){
  Coordinate3D agentsCoord = agent->getPosition();

  double newX,newY,newZ,diffX,diffY,diffZ;
  double oldX,oldY,oldZ,oldDiffX,oldDiffY,oldDiffZ;
  double virtualOldX,virtualOldY,virtualOldZ;
  bool lower,upper;

  Coordinate3D pointA = site->getLowerLimits();
  Coordinate3D pointB = site->getUpperLimits();

  switch(site->getNumberOfSpatialDimensions()){
    case 2:
      oldX = agentsCoord.x - moveVec->x;
      oldY = agentsCoord.y - moveVec->y;
      oldZ = 0.0;

      if(!((lower = (agentsCoord.x >= pointA.x)) && (upper = (agentsCoord.x <=pointB.x )))){

        if(!lower){
          diffX = pointA.x - agentsCoord.x;
          newX = pointA.x + diffX;

          oldDiffX = oldX - pointA.x;
          virtualOldX = pointA.x - oldDiffX;
        }else{
          diffX = pointB.x - agentsCoord.x ;
          newX = pointB.x + diffX;

          oldDiffX = pointB.x - oldX;
          virtualOldX = pointB.x + oldDiffX;
        }
      }else{
        newX = agentsCoord.x;
        virtualOldX = oldX;
      }

      if(!((lower = (agentsCoord.y>=pointA.y)) && (upper = (agentsCoord.y<=pointB.y)))){
        if(!lower){
          diffY = pointA.y - agentsCoord.y;
          newY = pointA.y + diffY;

          oldDiffY = oldY - pointA.y;
          virtualOldY = pointA.y - oldDiffY;
        }else{
          diffY = pointB.y - agentsCoord.y;
          newY = pointB.y + diffY;

          oldDiffY = pointB.y - oldY;
          virtualOldY = pointB.y + oldDiffY;
        }
      }else{
        newY = agentsCoord.y;
        virtualOldY = oldY;
      }

      newZ = agentsCoord.z;
      virtualOldZ = oldZ;

      break;


    case 3:
      oldX = agentsCoord.x - moveVec->x;
      oldY = agentsCoord.y - moveVec->y;
      oldZ = agentsCoord.z - moveVec->z;

      if(!((lower = (agentsCoord.x>=pointA.x)) && (upper = (agentsCoord.x<=pointB.x)))){

        if(!lower){
          diffX = pointA.x - agentsCoord.x;
          newX = pointA.x + diffX;

          oldDiffX = oldX - pointA.x;
          virtualOldX = pointA.x - oldDiffX;
        }else{
          diffX = pointB.x - agentsCoord.x ;
          newX = pointB.x + diffX;

          oldDiffX = pointB.x - oldX;
          virtualOldX = pointB.x + oldDiffX;
        }
      }else{
        newX = agentsCoord.x;
        virtualOldX = oldX;
      }

      if(!((lower = (agentsCoord.y>=pointA.y)) && (upper = (agentsCoord.y<=pointB.y)))){
        if(!lower){
          diffY = pointA.y - agentsCoord.y;
          newY = pointA.y + diffY;

          oldDiffY = oldY - pointA.y;
          virtualOldY = pointA.y - oldDiffY;
        }else{
          diffY = pointB.y - agentsCoord.y;
          newY = pointB.y + diffY;

          oldDiffY = pointB.y - oldY;
          virtualOldY = pointB.y + oldDiffY;
        }
      }else{
        newY = agentsCoord.y;
        virtualOldY = oldY;
      }

      if(!((lower = (agentsCoord.z>=pointA.z)) && (upper = (agentsCoord.z<=pointB.z)))){
        if(!lower){
          diffZ = pointA.z - agentsCoord.z;
          newZ = pointA.z + diffZ;

          oldDiffZ = oldZ - pointA.z;
          virtualOldZ = pointA.z - oldDiffZ;
        }else{
          diffZ =  pointB.z - agentsCoord.z;
          newZ = pointB.z + diffZ;

          oldDiffZ = pointB.z - oldZ;
          virtualOldZ = pointB.z + oldDiffZ;
        }
      }else{
        newZ = agentsCoord.z;
        virtualOldZ = oldZ;
      }
      break;


    default:
      oldX = agentsCoord.x - moveVec->x;
      oldY = agentsCoord.y - moveVec->y;
      oldZ = agentsCoord.z - moveVec->z;

      if(!((lower = (agentsCoord.x>=pointA.x)) && (upper = (agentsCoord.x<=pointB.x)))){

        if(!lower){
          diffX = pointA.x - agentsCoord.x;
          newX = pointA.x + diffX;

          oldDiffX = oldX - pointA.x;
          virtualOldX = pointA.x - oldDiffX;
        }else{
          diffX = pointB.x - agentsCoord.x ;
          newX = pointB.x + diffX;

          oldDiffX = pointB.x - oldX;
          virtualOldX = pointB.x + oldDiffX;
        }
      }else{
        newX = agentsCoord.x;
        virtualOldX = oldX;
      }

      if(!((lower = (agentsCoord.y>=pointA.y)) && (upper = (agentsCoord.y<=pointB.y)))){
        if(!lower){
          diffY = pointA.y - agentsCoord.y;
          newY = pointA.y + diffY;

          oldDiffY = oldY - pointA.y;
          virtualOldY = pointA.y - oldDiffY;
        }else{
          diffY = pointB.y - agentsCoord.y;
          newY = pointB.y + diffY;

          oldDiffY = pointB.y - oldY;
          virtualOldY = pointB.y + oldDiffY;
        }
      }else{
        newY = agentsCoord.y;
        virtualOldY = oldY;
      }

      if(!((lower = (agentsCoord.z>=pointA.z)) && (upper = (agentsCoord.z<=pointB.z)))){
        if(!lower){
          diffZ = pointA.z - agentsCoord.z;
          newZ = pointA.z + diffZ;

          oldDiffZ = oldZ - pointA.z;
          virtualOldZ = pointA.z - oldDiffZ;
        }else{
          diffZ =  pointB.z - agentsCoord.z;
          newZ = pointB.z + diffZ;

          oldDiffZ = pointB.z - oldZ;
          virtualOldZ = pointB.z + oldDiffZ;
        }
      }else{
        newZ = agentsCoord.z;
        virtualOldZ = oldZ;
      }
      break;

  }


  Coordinate3D prevMove{newX-virtualOldX,newY-virtualOldY,newZ-virtualOldZ};
  //corrected previous move in order to simulate the reflecting property of the border
  agent->getMovement()->setPreviousMove(&prevMove);
  Coordinate3D newPos = {newX,newY,newZ};
  auto newShift = newPos - agent->getPosition();
  //cout << "now reflecting boundary handling " << '\n';
  agent->shiftPosition(&newShift, current_time, agent->getSurface()->getBasicSphereOfThis(), "");
//  agent->setPosition(Coordinate3D{newX,newY,newZ});
}

std::string ReflectingBoundaries::getTypeName(){
  return "ReflectingBoundaries";
}

Coordinate3D ReflectingBoundaries::update_pos(Coordinate3D pos) {
    auto upper_limit = site->getUpperLimits();
    auto lower_limit = site->getLowerLimits();

    // Reflecting in x-direction
    if (pos.x < lower_limit.x){
        pos.x = lower_limit.x + (lower_limit.x - pos.x);
    }
    else if (pos.x > upper_limit.x){
        pos.x = upper_limit.x - (pos.x - upper_limit.x);
    }

    // Reflecting in y-direction
    if (pos.y < lower_limit.y){
        pos.y = lower_limit.y + (lower_limit.y - pos.y);
    }
    else if (pos.y > upper_limit.y){
        pos.y = upper_limit.y - (pos.y - upper_limit.y);
    }

    // Reflecting in z-direction
    if (pos.z < lower_limit.z){
        pos.z = lower_limit.z + (lower_limit.z - pos.z);
    }
    else if (pos.z > upper_limit.z){
        pos.z = upper_limit.z - (pos.z - upper_limit.z);
    }
    return pos;
}
