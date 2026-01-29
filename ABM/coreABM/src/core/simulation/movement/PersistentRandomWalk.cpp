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

#ifndef M_PI
#define M_PI    3.14159265358979323846f
#endif

#include "PersistentRandomWalk.h"
#include "core/simulation/Agent.h"
#include "core/simulation/Site.h"

PersistentRandomWalk::PersistentRandomWalk(Agent* agent, double persistenceTime , double speed, unsigned int spatial_dimensions) : Movement(spatial_dimensions) {
  this->agent_ = agent;
  this->persistence_time_ = persistenceTime;
  site_ = agent->getSite();

  sampler_ = std::make_unique<InversionSampler>();
  sampler_->setSampleFunction(0);

  if (agent->getInitialTime() > 0.0) {
    persistence_time_start_ = 0;
    persistence_time_left_ = site_->getRandomGenerator()->generateDouble() * persistenceTime;
  } else {
    persistence_time_start_ = site_->getRandomGenerator()->generateDouble() * persistenceTime;
    persistence_time_left_ = persistence_time_start_;
  }

  current_velocity_ = std::make_unique<Coordinate3D>();
  this->speed_ = speed;
  persistence_direction_ = std::make_unique<Coordinate3D>();
}


Coordinate3D* PersistentRandomWalk::move(double timestep, double diffusionConstant) {
  setCurrentTimestep(timestep);
  if (persistentMove()) {
    *current_velocity_ = *movePersistent(timestep);
  } else {
    moveRandomly(timestep);
    setNewPersistence();
  }
  decrementLeftTime(timestep);
  *current_move_ = *current_velocity_;
  return current_move_.get();
}

Coordinate3D* PersistentRandomWalk::movePersistent(double timestep){
  double length = persistence_direction_->getMagnitude();
  if (length == 0) {
    *persistence_direction_ = *moveRandomly(timestep);
  } else {
    *persistence_direction_ = site_->generatePersistentDirectionVector(*current_pos_,
                                                                       speed_ * timestep,
                                                                       *persistence_direction_,
                                                                       persistent_angle_alpha_2_d_);

  }

  return persistence_direction_.get();
}

Coordinate3D *PersistentRandomWalk::moveRandomly(double timestep) {

  double length = speed_ * timestep;
  *current_velocity_ = site_->generateRandomDirectionVector(*current_pos_, length);
  persistent_angle_alpha_2_d_ = site_->getLatestAlpha2dTurningAngle();

  return current_velocity_.get();
}

void PersistentRandomWalk::decrementLeftTime(double timestep){
  persistence_time_left_ -= timestep;
}

void PersistentRandomWalk::setNewPersistence() {
  *persistence_direction_ = *current_velocity_;
  persistence_time_left_ += persistence_time_;
}

bool PersistentRandomWalk::persistentMove() const {
  return (persistence_time_left_ > 0);
}

double PersistentRandomWalk::getStartingTime(){
  return round(persistence_time_start_);
}

void PersistentRandomWalk::setPreviousMove(Coordinate3D *prevMove) {
  *persistence_direction_ = *prevMove;
  persistent_angle_alpha_2_d_ = site_->getLatestAlpha2dTurningAngle();
}