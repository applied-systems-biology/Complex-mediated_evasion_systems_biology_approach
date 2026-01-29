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

#ifndef CORE_SIMULATION_INTERACTIONSTATE_H
#define CORE_SIMULATION_INTERACTIONSTATE_H

#include <utility>

#include "core/simulation/interactiontypes/InteractionType.h"
#include "core/simulation/Cell.h"

class Analyser;
class Interaction;

class InteractionState {
public:
  // Class for handling the interactions specified in the simulator-config
    InteractionState(std::string state_name, Interaction *interaction,
                     std::unique_ptr<InteractionType> interaction_type, bool end_state)
            : current_state_(std::move(state_name)), interaction_(interaction), end_state_(end_state),
              interaction_type_(std::move(interaction_type)) {}

    void addNextStateWithRate(const std::string &name_next_state, Rate *rate);
    void addNextStateWithRate(std::map<std::string, Rate *> next_states_rates);

    void fireInteractionEvent(const std::string &next_state, double current_time);
    void handleInteraction(Cell *cell, double timestep, double current_time);
    void stateTransition(double timestep, double current_time);
    [[nodiscard]] bool isEndState() const;
    [[nodiscard]] const std::string &getStateName() const;
    [[nodiscard]] std::string getInteractionType() const;

protected:
    void selectNextState(double timestep, Randomizer *randomizer, Condition *condition, Site *site, Cell *cell);
    bool end_state_{};
    std::string current_state_;
    std::string next_state_;
    Interaction *interaction_;
    std::map<std::string, Rate *> next_states_rates_;
    std::unique_ptr<InteractionType> interaction_type_;
};
#endif /* CORE_SIMULATION_INTERACTIONSTATE_H */
