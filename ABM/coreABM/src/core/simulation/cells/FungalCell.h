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


#ifndef CORE_SIMULATION_FUNGALCELL_H
#define CORE_SIMULATION_FUNGALCELL_H

#include "core/simulation/Cell.h"

class FungalCell : public Cell {
public:
    FungalCell(std::unique_ptr<Coordinate3D> c, int id, Site *site,
               double time_delta, double current_time) : Cell(std::move(c), id, site, time_delta, current_time) {}

    void doMorphologicalChanges(double timestep, double current_time);
    std::string getTypeName();
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters* parameters);
    abm::util::SimulationParameters::FungalParameters* getFungalParameter(){return fc_parameters;};

private:
    void handleInteractionEvent(InteractionEvent *ievent, double current_time) final;
    bool phagocytosed{};
    abm::util::SimulationParameters::FungalParameters* fc_parameters;
};

#endif /* CORE_SIMULATION_FUNGALCELL_H */
