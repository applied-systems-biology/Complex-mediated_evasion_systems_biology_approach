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

#ifndef DEFENSIVE_H
#define    DEFENSIVE_H

#include "core/simulation/Cell.h"

class Defensive : public Cell {
public:
  // Class for defensive cells
    Defensive(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Cell(std::move(c), id, site, time_delta, current_time) {};

    /*!
     * Sets up ImmuneCell from input parameters
     * @param time_delta Double that contains time step
     * @param current_time Double that contains current time
     * @param parameters SimulationParameters with agent parameters
     */
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) final;

    /// Returns "Defensive"
    std::string getTypeName() final;
    void update_nb_bounds();
    void setup_nb_bounds(int nb_bounds);
    int get_nb_bound(){return nb_bound_;}

private:
  void handleInteractionEvent(InteractionEvent *ievent, double current_time) final;
  int nb_bound_;
};

#endif    /* DEFENSIVE_H */

