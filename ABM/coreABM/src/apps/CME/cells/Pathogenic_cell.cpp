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

#include "Pathogenic_cell.h"
#include "core/simulation/Interactions.h"


std::string Pathogenic_cell::getTypeName() {
    return "Pathogenic_cell";
}

void Pathogenic_cell::setup(double time_delta,
                double current_time,
                abm::util::SimulationParameters::AgentParameters *parameters) {
    total_uptake = 0.0;

    positionShiftAllowed = false;
    auto param = static_cast<abm::utilCME::AgentsParametersCME*>(parameters);
    death_ = param->death_;
    for (auto& [name, data] : param->secretion_) {
        initial_rate_[name] = data.rate;
        //current_uptake_[data.mol_uptaken] = 0.0;
        //previous_uptake_[data.mol_uptaken].push_back(0.0);
        std::map<std::string, std::vector<Coordinate3D>> init;

        init[data.mol_uptaken];
        uptake_.push_back(init);
        last_time_uptake_updated[data.mol_uptaken] = 0.0;
        if (data.type != "constant"){
            data.rate = 0.0;
        }
        secretion_.emplace_back(std::pair<std::string, abm::utilCME::secretion_rate>(name, data));
    }
    //Set the fitted_proba_survival_relation_, form a text input file?? Or by "hand"??
    //fitted_proba_survival_relation_[100] = 8e-5;
    proba_survival_ = 1.0;
    Cell::setup(time_delta, current_time, parameters);
}

void Pathogenic_cell::update_total_uptake() {
    ++total_uptake;
}

void Pathogenic_cell::update_current_uptake(double current_time, double dt, std::string name, Coordinate3D pos){
    //DEBUG_STDOUT("UPTAKE HAPPENING: t= " << current_time);
    if (current_time > 0 && current_time == last_time_uptake_updated[name]){
        //current_uptake_[name] += 1;
        auto& current_uptake = uptake_.back();
        current_uptake[name].push_back(pos);
    }
    else{
        std::map<std::string, std::vector<Coordinate3D>> new_uptake;
        new_uptake[name].push_back(pos);
        uptake_.push_back(new_uptake);
        //current_uptake_[name] = 1;
        last_time_uptake_updated[name] = current_time;
    }
}

void Pathogenic_cell::move(double timestep, double current_time) {
}
void Pathogenic_cell::secreting_dynamics(double current_time, double dt) {
    for (const auto& [name, data]: secretion_) {
        if (data.type == "variable") {
            update_secretion_rate(current_time, dt);
        }
        if (data.rate != 0 && data.lag < current_time) {
            if (!data.spatial) {
                // SYSTEM_STDOUT("SECRETION RATE = " << data.rate * dt);
                dynamic_cast<CuboidSiteCME*>(site)->secretion_of_agents_at_cell_surface_randomly(dt, current_time, data.rate, name, *position, morphology_.get()->getBasicSphereOfThis()->getRadius());
            }
            else{
                auto& current_uptake = uptake_.front();
                dynamic_cast<CuboidSiteCME*>(site)->secretion_of_agents_at_cell_surface(dt, current_time, data.rate, name, *position, morphology_.get()->getBasicSphereOfThis()->getRadius(), current_uptake[data.mol_uptaken]);
            }
        }
    }
    //for (const auto& ag:site->getAgentManager()->getAllAgentTypes()){
    //    current_uptake_[ag] = 0.0;
    //}
}
void Pathogenic_cell::update_secretion_rate(double current_time, double dt) {
    for (auto& [name, data]: secretion_) {
        for (const auto& cell_uptaken: site->getAgentManager()->getAllAgentTypes()) {
            if (data.mol_uptaken == cell_uptaken) {
                //DEBUG_STDOUT("UPTAKE SIZE = " << uptake_.size());
                if (last_time_uptake_updated[cell_uptaken] < current_time) {
                    std::map<std::string, std::vector<Coordinate3D>> init;
                    init[data.mol_uptaken];
                    uptake_.push_back(init);
                }
                if (uptake_.size() > std::max(1.0, data.lag/dt)) { //last_time_uptake_updated[cell_uptaken]
                    //previous_uptake_[cell_uptaken].push_back(current_uptake_[cell_uptaken]);
                    //previous_uptake_[cell_uptaken].erase(previous_uptake_[cell_uptaken].begin());
                    uptake_.erase(uptake_.begin());
                }
                else {
                    //previous_uptake_[cell_uptaken].push_back(0.0);

                    std::map<std::string, std::vector<Coordinate3D>> init;
                    init[data.mol_uptaken];
                    uptake_.push_back(init);
                }
                if (!uptake_.empty()) {
                    data.rate = initial_rate_[name] * uptake_.front()[cell_uptaken].size() / dt;
                }
                else{
                    data.rate = 0.0;
                }
            }
        }
    }
}
abm::utilCME::secretion_rate& Pathogenic_cell::get_secretion(std::string name) {
    for (auto& sec: secretion_){
        if (sec.first == name){
            return sec.second;
        }
    }
    throw std::runtime_error("Name not found in secretion vector. Please check your screening param name in the config.json file");
}

void Pathogenic_cell::doAllActionsForTimestep(double timestep, double current_time) {
    secreting_dynamics(current_time, timestep);
}