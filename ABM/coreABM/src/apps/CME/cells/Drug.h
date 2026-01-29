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

#ifndef COREABM_DRUG_H
#define COREABM_DRUG_H

#include "core/simulation/Cell.h"

class Drug : public Cell {
  public:
    // Class for Drug agents
    Drug(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Cell(std::move(c), id, site, time_delta, current_time) {};

    /*!
     * Sets up Drug from input parameters
     * @param time_delta Double that contains time step
     * @param current_time Double that contains current time
     * @param parameters SimulationParameters with agent parameters
     */
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) final;

    /// Returns "Drug"
    std::string getTypeName() final;

    void update_nb_bounds();
    void setup_nb_bounds(int nb_bounds);
    int get_nb_bound(){return nb_bound_;}
    int get_nb_receptors_blocked(){return nb_receptors_blocked_;}

  private:
    void handleInteractionEvent(InteractionEvent *ievent, double current_time) final;
    int nb_bound_;
    int nb_receptors_blocked_;
};
#endif // COREABM_DRUG_H
