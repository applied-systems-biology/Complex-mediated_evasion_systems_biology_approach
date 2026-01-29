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

#include "AperiodicBoundaries.h"

#include "core/analyser/Analyser.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/AgentManager.h"
#include "core/simulation/Cell.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Interactions.h"
#include "core/simulation/Site.h"
#include "core/simulation/cells/ImmuneCell.h"
#include "core/simulation/states/InteractionState.h"

AperiodicBoundaries::AperiodicBoundaries(Site* site) : BoundaryCondition(site) {
}

void AperiodicBoundaries::handleBoundaryCross(Agent* agent, Coordinate3D* moveVec, double current_time){
  std::string agentType = agent->getTypeName();

  //check for interacting agents e.g. ingestion/ phagocyte-fungus-Interaction

  const auto interactions = agent->getInteractions()->getAllInteractions();
  bool agentIsPartOfComplex = false;
  std::vector<Agent*> agentComplex = agent->getComplex();
  if(agentComplex.size() > 1){
    agentIsPartOfComplex = true;
  }

  if(agentIsPartOfComplex){
    site->getAgentManager()->insertExistingComplexAtBoundary(site, agentComplex, agent->getId(), current_time);
  }else{

    if(abm::util::isSubstring("FungalCell", agentType)){
      std::string stateName = agent->getCurrentCellState()->getStateName();
      site->getAgentManager()->replaceAgent(site, agent, 0, current_time);
      site->getAgentManager()->insertExistingAgentAtBoundary(site, agentType, stateName, current_time);

    }else{
      site->getAgentManager()->replaceAgent(site, agent, 0, current_time);
      //auto has_contact = dynamic_cast<ImmuneCell*>(agent)->has_made_contact();
      site->getAgentManager()->insertAgentAtBoundary(site, agentType, current_time);
    }
  }
//
//    default:
//      site->getAgentManager()->replaceAgent(site, agent, 0, current_time);
//      site->getAgentManager()->insertAgentAtBoundary(site, agentType, current_time);
//      break;
}


std::string AperiodicBoundaries::getTypeName(){
  return "AperiodicBoundaries";
}

Coordinate3D AperiodicBoundaries::aperiodicMovementTransformation(Coordinate3D* prevMove, Coordinate3D* boundaryPoint, double delta_time){
  double prevSpeed = prevMove->getMagnitude();

  Coordinate3D virtualPoint = Coordinate3D();
  RandomWalk rw = RandomWalk(site,prevSpeed/delta_time);
  rw.setSite(site);
  rw.setCurrentPosition(boundaryPoint);

  Coordinate3D newPrevMove = Coordinate3D();

  do{
    virtualPoint = *boundaryPoint;
    newPrevMove = *rw.calculateRandomMove(delta_time);
    virtualPoint += newPrevMove;
  }while(!site->containsPosition(virtualPoint));

  return newPrevMove;
}
Coordinate3D AperiodicBoundaries::update_pos(Coordinate3D pos) {
    return site->getRandomBoundaryPoint();
}
