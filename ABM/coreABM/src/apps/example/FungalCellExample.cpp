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


#include "FungalCellExample.h"
#include "core/simulation/Site.h"
#include "core/simulation/morphology/SphericalMorphology.h"

#include "core/utils/macros.h"

void FungalCellExample::doMorphologicalChanges(double timestep, double current_time) {
}

void FungalCellExample::setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters * parameters) {
    fc_parameters = static_cast<abm::utilExample::FungalParametersExample * >(parameters);
    std::cout << "Set example parameter for FungalCellExample id=" << getId() << ": hyphal branch: " << fc_parameters->hyphal_growth << " - ";
    std::cout << "example: " << fc_parameters->other_example_parameter << "\n";
    Cell::setup(time_delta, current_time, parameters);
}
