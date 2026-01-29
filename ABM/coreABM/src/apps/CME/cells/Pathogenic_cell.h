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

//
// Created by ybachelot on 20.04.23.
//

#ifndef COREABM_PATHOGENIC_CELL_H
#define COREABM_PATHOGENIC_CELL_H

#include "apps/CME/CuboidSiteCME.h"
#include "apps/CME/io_utilsCME.h"
#include "core/simulation/Cell.h"

class Pathogenic_cell : public Cell{
  public:

    Pathogenic_cell(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Cell(std::move(c), id, site, time_delta, current_time) {
        this->site = dynamic_cast<CuboidSiteCME*>(site);
          };

    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) final;

    void update_total_uptake();

    int get_total_uptake(){return total_uptake;}

    std::string getTypeName() final;

    void move(double timestep, double current_time) final;

    void secreting_dynamics(double current_time, double dt);

    void update_secretion_rate(double current_time, double dt);

    void update_current_uptake(double current_time, double dt, std::string name, Coordinate3D pos);

    bool get_death(){return death_;}
    abm::utilCME::secretion_rate& get_secretion(std::string name);

    double get_proba_survival(){return proba_survival_;}

    void doAllActionsForTimestep(double timestep, double current_time) override;


protected:
    int total_uptake{};
    bool death_{};
    std::vector<std::pair<std::string, abm::utilCME::secretion_rate>> secretion_;
    //std::map<std::string, double> current_uptake_;
    //std::map<std::string, std::vector<double>> previous_uptake_;
    std::map<std::string, double> initial_rate_;
    std::map<std::string, double> last_time_uptake_updated;
    std::vector<std::map<std::string, std::vector<Coordinate3D>>> uptake_;
    std::map<double, double> fitted_proba_survival_relation_;
    double proba_survival_;
};

#endif // COREABM_PATHOGENIC_CELL_H
