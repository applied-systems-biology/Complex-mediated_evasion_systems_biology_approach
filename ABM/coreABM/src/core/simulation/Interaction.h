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

#ifndef CORE_SIMULATION_INTERACTION_H
#define CORE_SIMULATION_INTERACTION_H

#include <memory>
#include <map>
#include <queue>

#include "core/simulation/Condition.h"
#include "core/simulation/cells/interaction/InteractionEvent.h"
#include "core/simulation/neighbourhood/Collision.h"
#include "core/analyser/Analyser.h"

class Cell;
class InteractionState;

class Interaction {
public:
  // Class for providing the key functionality of an interaction.
  Interaction(std::string identifier, Cell *cell1, Cell *cell2, double time_delta, double current_time);
  virtual ~Interaction() = default;

  virtual void handle(Cell *cell, double timestep, double current_time);
  [[nodiscard]] virtual std::string getInteractionName() const;
  std::shared_ptr<Collision> getNextCollision();
  Condition *getCurrentCondition();
  InteractionState *getCurrentState() { return interactionState.get(); };
  Cell *getFirstCell();
  Cell *getSecondCell();
  Cell *getOtherCell(Cell *cell);
  void setDelted() { setDelete = true; };
  void setInitialState(double time_delta, double current_time, Cell *cell = nullptr);
  void setState(std::string nameOfState);
  void fireInteractionEvent(InteractionEvent *ievent, double current_time);
  void addCurrentCollision(std::shared_ptr<Collision> collision) {
    this->currentCollisions.push(std::move(collision));
  };
  bool isActive();
  bool isDelted() const { return setDelete; };
  std::string getIdentifier(){return identifier_;}
  void close();

protected:
  unsigned int interactionId;
  Cell *cellOne;
  Cell *cellTwo;
  std::map<Cell *, std::shared_ptr<Condition>> cellularConditions;
  std::map<std::string, std::shared_ptr<Condition>> cellularStringConditions;
  std::queue<std::shared_ptr<Collision>> currentCollisions;
  std::shared_ptr<InteractionState> interactionState;
  std::shared_ptr<InteractionState> oldinteractionState;
  Condition *currentCondition;

private:
  bool isActiven;
  bool setDelete;
  std::string identifier_;
};

#endif /* CORE_SIMULATION_INTERACTION_H */
