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

#include "Drug.h"
#include "apps/CME/io_utilsCME.h"

std::string Drug::getTypeName() {
    return "Drug";
}

void Drug::setup(double time_delta,
                double current_time,
                abm::util::SimulationParameters::AgentParameters *parameters) {

    auto *param = static_cast<abm::utilCME::Drug *>(parameters);
    nb_receptors_blocked_ = param->nb_recep_blocked;
    Cell::setup(time_delta, current_time, parameters);
}

void Drug::handleInteractionEvent(InteractionEvent *ievent, double current_time) {
    if (ievent->getNextState() == "Complex") {
        update_nb_bounds();
        setDeleted();
    }
}
void Drug::update_nb_bounds() {
    ++nb_bound_;
}
void Drug::setup_nb_bounds(int nb_bounds) {
    nb_bound_ = nb_bounds;
}