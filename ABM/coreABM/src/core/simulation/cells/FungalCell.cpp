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


#include "FungalCell.h"
#include "core/simulation/Site.h"
#include "core/simulation/morphology/SphericalMorphology.h"

#include "core/utils/macros.h"


std::string FungalCell::getTypeName() {
  return "FungalCell";
}

void FungalCell::handleInteractionEvent(InteractionEvent *ievent, double current_time) {

    if (ievent->getNextState() == "Phagocytose") {
        auto fState = getCellStateByName("FungalPhagocytosed");
        if (fState != 0) {
            setState(fState);
            INFO_STDOUT("Phagocytose: FungalCell was phagocytosed");
            phagocytosed = true;
            getSite()->getAgentManager()->removeFungalCellFromList(this->getId(), current_time);
        }
    }
    if (ievent->getNextState() == "Lysis") {
        auto fState = getCellStateByName("Death");
        if (fState != 0) {
            setState(fState);
            INFO_STDOUT("Lysis: Death of FungalCell");
        }
        setDeleted();
    }
}

void FungalCell::doMorphologicalChanges(double timestep, double current_time) {
}

void FungalCell::setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters * parameters) {
    fc_parameters = static_cast<abm::util::SimulationParameters::FungalParameters * >(parameters);
    Cell::setup(time_delta, current_time, parameters);
}
