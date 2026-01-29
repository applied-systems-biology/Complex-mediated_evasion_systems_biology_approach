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

#ifndef CORE_SIMULATION_IMMUNECELL_H
#define CORE_SIMULATION_IMMUNECELL_H

#include "core/simulation/Cell.h"

class ImmuneCell : public Cell {
public:
  // Class for immune cells
    ImmuneCell(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Cell(std::move(c), id, site, time_delta, current_time) {};

    /*!
     * Sets up ImmuneCell from input parameters
     * @param time_delta Double that contains time step
     * @param current_time Double that contains current time
     * @param parameters SimulationParameters with agent parameters
     */
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters);

    /// Returns "ImmuneCell"
    std::string getTypeName();
    void move(double timestep, double current_time) final;

private:
    void handleInteractionEvent(InteractionEvent *ievent, double current_time);
    unsigned int currentNoOfUptakes;

    double radius;
    Coordinate3D cumulativePersistenceGradient;
    abm::util::SimulationParameters::ImmuneCellParameters* ic_parameters;
};

#endif /* CORE_SIMULATION_IMMUNECELL_H */
