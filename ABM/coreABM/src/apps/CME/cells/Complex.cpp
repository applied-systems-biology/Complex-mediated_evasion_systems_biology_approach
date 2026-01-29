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

#include "Complex.h"

#include "apps/CME/cells/Drug.h"
#include "apps/CME/io_utilsCME.h"
#include "core/analyser/Analyser.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/Site.h"
#include "core/simulation/factories/CellFactory.h"
#include "core/simulation/Interactions.h"

#include <cmath>

std::string Complex::getTypeName() {
    return "Complex";
}



void Complex::setup(double time_delta,
                       double current_time,
                       abm::util::SimulationParameters::AgentParameters *parameters) {

    lifetime_ = 0.0;
    auto *param = static_cast<abm::utilCME::Complex *>(parameters);
    total_nb_receptors_ = param->nb_receptors;
    free_receptors_ = total_nb_receptors_;
    Cell::setup(time_delta, current_time, parameters);
}

void Complex::update_lifetime(){
    ++lifetime_;
}

void Complex::unbinding_dynamics(double timestep, double current_time, std::string mol_to_remove, int nb_recep_freed) {
    auto possible = false;
    auto it = std::find(molecules_in_complex_.begin(), molecules_in_complex_.end(), mol_to_remove);
    if (it != molecules_in_complex_.end()) {
        possible = true;
    }
    if (possible){
        //DEBUG_STDOUT("Unbinding EVENT: " << mol_to_remove);

        auto *agent_manager = site->getAgentManager();
        auto initial_pos = *position;
        auto radius_comp = this->getMorphology()->getBasicSphereOfThis()->getRadius();
        bool found_r1 = false;
        bool found_r2 = false;
        if (molecules_in_complex_.size() == 0){
            ERROR_STDERR("This should not happen");
        }
        else if (molecules_in_complex_.size() <= 1) {
            auto radius_1 = 0.00104, radius_2 = 0.00334;
            for (const auto& ag: agent_manager->getAllAgents()){
                if (ag->getTypeName() == mol_to_remove){
                    radius_1 = ag->getMorphology()->getBasicSphereOfThis()->getRadius();
                    found_r1= true;
                }
                else if (ag->getTypeName() == "Defensive"){
                    radius_2 = ag->getMorphology()->getBasicSphereOfThis()->getRadius();
                    found_r2 = true;
                }
                if (found_r1 && found_r2){
                    break;
                }
            }
            auto min_sep = radius_1 + radius_2;
            auto max_sep = radius_1 + radius_2 + std::min(radius_1, radius_2);

            Coordinate3D pos1, pos2;

            auto theta = site->getRandomGenerator()->generateDouble(0, M_PI);
            auto phi1 = site->getRandomGenerator()->generateDouble(0, 2 * M_PI);
            auto phi2 = site->getRandomGenerator()->generateDouble(0, 2 * M_PI);

            // Select a random distance between min_sep and max_sep
            auto r = site->getRandomGenerator()->generateDouble(min_sep, max_sep);

            pos1.x = initial_pos.x + r * 0.5 * sin(theta) * cos(phi1);
            pos1.y = initial_pos.y + r * 0.5 * sin(theta) * sin(phi1);
            pos1.z = initial_pos.z + r * 0.5 * cos(theta);

            pos2.x = initial_pos.x - r * 0.5 * sin(theta) * cos(phi2);
            pos2.y = initial_pos.y - r * 0.5 * sin(theta) * sin(phi2);
            pos2.z = initial_pos.z - r * 0.5 * cos(theta);


            if (!site->containsPosition(pos1)){
                pos1 = site->get_boundary()->update_pos(pos1);
            }
            if (!site->containsPosition(pos2)){
                pos2 = site->get_boundary()->update_pos(pos2);
            }
            agent_manager->emplace_back(site->getCellFactory()->createCell(mol_to_remove,
                                                                           std::make_unique<Coordinate3D>(pos1),
                                                                           agent_manager->generateNewID(), site, timestep, current_time, ""));
            agent_manager->emplace_back(site->getCellFactory()->createCell("Defensive",
                                                                           std::make_unique<Coordinate3D>(pos2),
                                                                           agent_manager->generateNewID(), site, timestep, current_time, ""));
            site->getMeasurments()->addValues<HistogramMeasurement>("lifetime", std::make_pair("timestep", timestep), std::make_pair("lifetime_complex", lifetime_));
            remove_molecule_from_complex(mol_to_remove);
            setDeleted();
        }
        else{
            auto radius_1 = 0.0;
            for (const auto& ag: agent_manager->getAllAgents()){
                if (ag->getTypeName() == mol_to_remove) {
                    radius_1 = ag->getMorphology()->getBasicSphereOfThis()->getRadius();
                    break;
                }
            }
            auto min_sep = radius_1 + radius_comp;
            auto max_sep = 2*radius_1 + radius_comp;

            Coordinate3D pos1;

            auto position_valid = false;

            do {
                // Generate a random distance within the specified range
                double r = site->getRandomGenerator()->generateDouble(min_sep, max_sep);

                // Generate random angles for the spherical coordinates
                double theta = site->getRandomGenerator()->generateDouble(0, M_PI);
                double phi = site->getRandomGenerator()->generateDouble(0, 2 * M_PI);

                // Convert spherical coordinates to Cartesian coordinates to get the child's position
                pos1.x = initial_pos.x + r * sin(theta) * cos(phi);
                pos1.y = initial_pos.y + r * sin(theta) * sin(phi);
                pos1.z = initial_pos.z + r * cos(theta);

                // Check if the position is within the bounds of the simulation space
                position_valid = site->containsPosition(pos1);

            } while (!position_valid);
            remove_molecule_from_complex(mol_to_remove);
            increase_free_receptor(nb_recep_freed);
            agent_manager->emplace_back(site->getCellFactory()->createCell(mol_to_remove,
                                                                           std::make_unique<Coordinate3D>(pos1),
                                                                           agent_manager->generateNewID(), site, timestep, current_time, ""));

            setExistingState("InitialCellState", timestep, current_time);
        }
    }
    else {
        setExistingState("InitialCellState", timestep, current_time);
    }
    update_radius_diffusion();
}

void Complex::remove_molecule_from_complex(std::string mol) {
    auto it = std::find(molecules_in_complex_.begin(),
                        molecules_in_complex_.end(),
                        mol);
    if (it != molecules_in_complex_.end())
        molecules_in_complex_.erase(it);
}

void Complex::add_molecule_to_complex(std::string mol) {
    if (mol != ""){
        molecules_in_complex_.push_back(mol);
    }
}
void Complex::increase_free_receptor(int nb) {
    free_receptors_ += nb;
    if (free_receptors_ > total_nb_receptors_){
        free_receptors_ = total_nb_receptors_;
    }
}
void Complex::decrease_free_receptor(int nb) {
    free_receptors_ -= nb;
    if (free_receptors_ < 0){
        free_receptors_ = 0;
    }
}
int Complex::get_nb_of_mol_in_complex(std::string mol) {
    auto res = 0;
    for (auto &name: molecules_in_complex_){
        if (name == mol){
            ++res;
        }
    }
    return res;
}
void Complex::extension(Cell& cell, int nb_recep_block) {
    //update free_receptors, content_complex, ...
    decrease_free_receptor(nb_recep_block);
    add_molecule_to_complex(cell.getTypeName());
    update_radius_diffusion();
}

void Complex::doAllActionsForTimestep(double timestep, double current_time) {
    enum class Task {
        MOVEMENT,
        INTERACTIONS_AND_STATES,
        MORPHOLOGY_CHANGE,
        MOLECULE_INTERACTION
    };
    enum class Change {
        STATES,
        INTERACTIONS
    };

    // Only treat once per timestep
    if (agentTreatedInCurrentTimestep(current_time)) {
        if (!is_deleted_) move(timestep, current_time);
    } else {

        // Generate random order of tasks
        std::vector<unsigned> perm = abm::util::generateRandomPermutation(
                site->getRandomGenerator(), 4
        );

        for (auto idx: perm) {
            switch (static_cast<Task>(idx)) {
                case Task::MOVEMENT:
                    if (!isDeleted()) {
                        move(timestep, current_time);
                    }
                    break;

                case Task::INTERACTIONS_AND_STATES: {
                    // Randomize STATES vs INTERACTIONS
                    std::vector<unsigned> statePerm = abm::util::generateRandomPermutation(
                            site->getRandomGenerator(), 2
                    );
                    for (auto cidx: statePerm) {
                        switch (static_cast<Change>(cidx)) {
                            case Change::STATES:
                                if (!isDeleted()) {
                                    cellState->stateTransition(timestep, current_time);
                                    if (cellState->checkForDeath(current_time)) {
                                        setDeleted();
                                        return;
                                    }
                                }
                                break;

                            case Change::INTERACTIONS:
                                if (!isDeleted()) {
                                    interactions->doWholeProcess(
                                            timestep, current_time, site->getMeasurments()
                                    );
                                }
                                break;
                        }
                    }
                }
                    break;

                case Task::MORPHOLOGY_CHANGE:
                    if (!isDeleted()) {
                        doMorphologicalChanges(timestep, current_time);
                    }
                    break;

                case Task::MOLECULE_INTERACTION:
                    if (!isDeleted()) {
                        // Complex-specific unbinding logic
                        const std::string stateName = getCurrentCellState()->getStateName();
                        if (stateName == "Unbound_AMP") {
                            unbinding_dynamics(timestep, current_time, "AMP", 1);
                        } else if (stateName == "Unbound_Drug") {
                            int nb_recep_freed = 0;
                            for (auto &agt: site->getAgentManager()->getAllAgents()) {
                                if (agt->getTypeName() == "Drug") {
                                    if (auto drug = dynamic_cast<Drug *>(agt.get())) {
                                        nb_recep_freed = drug->get_nb_receptors_blocked();
                                    }
                                    break;
                                }
                            }
                            unbinding_dynamics(timestep, current_time, "Drug", nb_recep_freed);
                        }
                    }
                    break;
            }
            if (isDeleted()) {
                break;
            }
        }

        // After all tasks, mark as treated
        setTimestepLastTreatment(current_time);

        // Cleanup spatial spheres if deleted
        if (isDeleted()) {
            for (const auto &sphere: getMorphology()->getAllSpheresOfThis()) {
                site->getNeighbourhoodLocator()->removeSphereRepresentation(sphere);
                site->getAgentManager()->removeSphereRepresentation(sphere);
            }
        }
    }
}

void Complex::update_radius_diffusion() {
    // update radius
    int count_AMP = std::count(molecules_in_complex_.begin(), molecules_in_complex_.end(), "AMP");
    int count_drug = std::count(molecules_in_complex_.begin(), molecules_in_complex_.end(), "Drug");
    auto mw_AMP = 4e3;

    auto drug_size = 0;
    if (count_drug > 0) {
        drug_size = (free_receptors_ - count_AMP) / count_drug;
    }

    auto mw_drug = drug_size * 4e3; // dependent on nb_receptor_blocked by Drug.

    auto mw_def = 130e3;
    auto molecular_weight = count_AMP * mw_AMP + count_drug * mw_drug + mw_def;

    auto radius = 0.066 * std::pow(molecular_weight, (1.0/3.0))*1e-3;
    //SYSTEM_STDOUT("radius = " << radius);

    morphology_->getBasicSphereOfThis()->setRadius(radius);

    // update diffusion
    auto eta = 4.5e-3;
    auto kb = 1.380649e-23;
    auto T = 310;
    auto D = (kb*T)/(radius*1e-6*6*M_PI*eta)*1e12;
    auto speed = sqrt(6*D/getSite()->getTimeStepping());
    movement_->setSpeed(speed);
}

