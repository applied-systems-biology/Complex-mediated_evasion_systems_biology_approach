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

#ifndef COMPLEX_H
#define    COMPLEX_H

#include "core/simulation/Cell.h"

class Complex : public Cell {
public:
  // Class for complex cells
    Complex(std::unique_ptr<Coordinate3D> c, int id, Site *site, double time_delta, double current_time)
        : Cell(std::move(c), id, site, time_delta, current_time) {};

    /*!
     * Sets up Complex from input parameters
     * @param time_delta Double that contains time step
     * @param current_time Double that contains current time
     * @param parameters SimulationParameters with agent parameters
     */
    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters *parameters) final;

    /// Returns "Complex"
    std::string getTypeName() final;

    void update_lifetime();
    void decrease_free_receptor(int nb);
    void increase_free_receptor(int nb);

    int get_total_nb_receptors() const {return total_nb_receptors_;};

    void unbinding_dynamics(double timestep, double current_time, std::string mol_to_remove, int nb_recep_freed);
    double get_lifetime() const {return lifetime_;}

    int get_free_receptors(){return free_receptors_;}
    void add_molecule_to_complex(std::string mol);
    void remove_molecule_from_complex(std::string mol);

    std::vector<std::string> get_complex_content(){return molecules_in_complex_;};

    int get_nb_of_mol_in_complex(std::string mol);
    void set_free_receptors(int nb){free_receptors_ = nb;};

    void extension(Cell& cell, int nb_recp_blocked);

    void doAllActionsForTimestep(double timestep, double current_time) override;

    void update_radius_diffusion();


private:
    double lifetime_;
    int total_nb_receptors_;

    std::vector<std::string> molecules_in_complex_;
    int free_receptors_;
};

#endif    /* COMPLEX_H */

