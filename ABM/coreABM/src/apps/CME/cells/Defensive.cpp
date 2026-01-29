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

#include "Defensive.h"

#include "core/analyser/Analyser.h"
#include "core/simulation/Site.h"


std::string Defensive::getTypeName() {
    return "Defensive";
}

void Defensive::setup(double time_delta,
                       double current_time,
                       abm::util::SimulationParameters::AgentParameters *parameters) {
    Cell::setup(time_delta, current_time, parameters);
}

void Defensive::handleInteractionEvent(InteractionEvent *ievent, double current_time){
    if (ievent->getNextState() == "Complex"){
        update_nb_bounds();
        setDeleted();
    }
}
void Defensive::update_nb_bounds() {
    ++nb_bound_;
}
void Defensive::setup_nb_bounds(int nb_bounds) {
    nb_bound_ = nb_bounds;
}
