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

#include <numeric>
#include <algorithm>

#include "core/simulation/Site.h"
#include "core/simulation/boundary-condition/AbsorbingBoundaries.h"
#include "core/simulation/boundary-condition/ReflectingBoundaries.h"
#include "core/simulation/boundary-condition/AperiodicBoundaries.h"
#include "core/simulation/boundary-condition/PeriodicBoundaries.h"
#include "core/simulation/AgentManager.h"
#include "core/utils/macros.h"
#include "external/json.hpp"
#include "core/utils/macros.h"

using json = nlohmann::json;

Site::Site(Randomizer *random_generator,
           std::shared_ptr<InSituMeasurements> measurements,
           std::string config_path, std::unordered_map<std::string, std::string> cmd_input_args,
           std::string output_path) : random_generator_(
        random_generator), measurements_(std::move(measurements)) {

    measurements_->setSite(this);
}

void Site::setBoundaryCondition() {
    if (boundary_type_ == "ReflectingBoundaries") {
        boundary_condition_ = std::make_unique<ReflectingBoundaries>(this);
    } else if (boundary_type_ == "AbsorbingBoundaries") {
        boundary_condition_ = std::make_unique<AbsorbingBoundaries>(this);
    } else if (boundary_type_ == "AperiodicBoundaries") {
        boundary_condition_ = std::make_unique<AperiodicBoundaries>(this);
    } else if (boundary_type_ == "PeriodicBoundaries") {
        boundary_condition_ = std::make_unique<PeriodicBoundaries>(this);
    } else {

        ERROR_STDERR("Boundary condition not found. Edit your simulator-config.json and use an existing boundary condition.");
        exit(1);
    }
}

void Site::doAgentDynamics(Randomizer *random_generator, SimulationTime &time) {

    const auto current_time = time.getCurrentTime();
    const auto dt = time.getCurrentDeltaT();

    const auto &all_agents = agent_manager_->getAllAgents();
    if (!all_agents.empty()) {
        // Loop over all agents (random order)
        const auto current_order = abm::util::generateRandomPermutation(random_generator, all_agents.size());
        for (auto agent_idx = current_order.begin(); agent_idx < current_order.end(); ++agent_idx) {
            auto curr_agent = all_agents[*agent_idx];
            if (nullptr != curr_agent) {
                // Do all actions for one timestep for each agent (-> Cell.cpp)
                curr_agent->doAllActionsForTimestep(dt, current_time);
                // Remove spherical representations if the current agent got deleted
                if (curr_agent->isDeleted()) {
                    for (const auto &sphere: curr_agent->getMorphology()->getAllSpheresOfThis()) {
                        neighbourhood_locator_->removeSphereRepresentation(sphere);
                        agent_manager_->removeSphereRepresentation(sphere);
                    }
                    curr_agent = nullptr;
                }
            }
        }

        // Clean up agents
        agent_manager_->cleanUpAgents(current_time);
        agent_manager_->inputOfAgents(current_time, random_generator);
        measurements_->observeMeasurements(time);
    }
}

bool Site::checkForStopping() const {
    return false;//stopSimulation;
}

//void Site::setStoppingCondition(const std::vector<std::string> &stopping_crit) {
//    for (auto criterion: stopping_crit) {
//        stopping_criteria.emplace_back(criterion);
//    }
//}

//void Site::stopRunForCertainState(Cell &cell, std::string state, double current_time) {
//    // Ends a simulation if all fungi were touched at least once (FTP)
//    // Cumulated first passage time (FTP) is clearance time (CT)
//    if (abm::util::isSubstring("FungalCell", cell.getTypeName())) {
//        if (state == "FungalPhagocytosed" or state == "Death") {
//            int cellid = cell.getId();
//            int number_of_previously_found_fungal_cells = detected_fungi_id.size();
//            detected_fungi_id.emplace_back(cellid);
//            std::sort(detected_fungi_id.begin(), detected_fungi_id.end());
//            detected_fungi_id.erase(std::unique(detected_fungi_id.begin(), detected_fungi_id.end()), detected_fungi_id.end());
//            if (number_of_previously_found_fungal_cells < detected_fungi_id.size()) {
//                DEBUG_STDOUT("FungalCell with id=" << cellid << " was detected by ImmuneCell at " << current_time);
//                int fungal_cells_remaining = agent_manager_->getInitFungalQuantity() - detected_fungi_id.size();
//                bool stopping_FPT = find(stopping_criteria.begin(), stopping_criteria.end(), "FirstPassageTime") != stopping_criteria.end();
//                if (stopping_FPT && fungal_cells_remaining == 0 && agent_manager_->getInitFungalQuantity() > 0) {
//                    stopSimulation = true;
//                    DEBUG_STDOUT("All FungalCells were detected.. end simulation");
//                }
//            }
//        }
//    }
//}

void Site::handleCmdInputArgs(std::unordered_map<std::string, std::string>cmd_input_args) {
    // Parameters to be screened or from cmd input
    // You can add your parameters you wanna screen below here
    if (cmd_input_args.size() > 0) {
        for (const auto &agent: parameters_.site_parameters->agent_manager_parameters.agents) {
            if (abm::util::isSubstring("ImmuneCell", agent->type)) {
                for (const auto &[key, value]: cmd_input_args) {
                    if ("icNum" == key) {
                        agent->number = std::stod(value);
                    }
                }
            }
            if (abm::util::isSubstring("FungalCell", agent->type)) {
                auto *fc_parameters = static_cast<abm::util::SimulationParameters::AgentParameters *>(agent.get());
                for (const auto &[key, value]: cmd_input_args) {
                    if ("fcNum" == key) {
                        fc_parameters->number = std::stoi(value);
                    }
                }
            }
        }
    }
}

void Site::receiveFrontendParameter(abm::util::SimulationParameters &sim_para, abm::util::InputParameters &inp_para,
                                    const std::string &config_path, const std::string &output_path, std::string sid) {

    const auto visualizer_config = static_cast<boost::filesystem::path>(config_path).append("visualisation-config.json").string();
    const auto fe_config_sid = static_cast<boost::filesystem::path>(output_path).append("frontend-api"+sid+".json").string();

    if (boost::filesystem::exists(fe_config_sid) && sid != "") {
        std::ifstream json_file_fe(fe_config_sid);
        json fe_para;
        json_file_fe >> fe_para;
        auto visp = abm::util::getViualizerParameters(visualizer_config);

        sim_para.max_time = fe_para["frontend-api"]["simulation_parameter"]["max_time"].value("value", 100);
        sim_para.time_stepping = fe_para["frontend-api"]["simulation_parameter"]["time_step"].value("value", 0.1);
//        auto csite = abm::util::SimulationParameters::CuboidSiteParameters{*sim_para.site_parameters};
        auto number_images = fe_para["frontend-api"]["simulation_parameter"]["output_images"].value("value", 10);

        if (number_images > 0) {
            visp.output_interval = int(sim_para.max_time / (number_images * sim_para.time_stepping));
        } else {
            visp.pov_active = false;
        }

        visp.output_interval = std::max(int(sim_para.max_time / (number_images * sim_para.time_stepping)), 1);
        visp.output_video = false; // important for the frontend, this must be false
        sim_para.dimensions = int(fe_para["frontend-api"]["simulation_parameter"]["dimensions"].value("value", 2));

        if (sim_para.dimensions == 3) {
            visp.camera_position = Coordinate3D{-visp.camera_position.z, visp.camera_position.z/4, visp.camera_position.z/2};
            if (sim_para.site_parameters->type == "CuboidSite") {
                auto* cuboid_parameters = static_cast<abm::util::SimulationParameters::CuboidSiteParameters*>(sim_para.site_parameters.get());
                cuboid_parameters->upper_bound.z = (cuboid_parameters->upper_bound.x + cuboid_parameters->upper_bound.y)/2;
                cuboid_parameters->lower_bound.z = (cuboid_parameters->lower_bound.x + cuboid_parameters->lower_bound.y)/2;
            }
        }

        for (auto& agent: sim_para.site_parameters->agent_manager_parameters.agents){
            if (agent->type == "FungalCell") {
                agent->number = int(fe_para["frontend-api"]["fungal_cell"]["number"].value("value", 10));
                agent->morphology_parameters.radius = fe_para["frontend-api"]["fungal_cell"]["radius"].value("value", 1.25);
                bool diff_colors = fe_para["frontend-api"]["fungal_cell"]["use_different_colors"].value("value", false);
                agent->morphology_parameters.color = diff_colors ? "random" : "redTransp";
            } else if (agent->type == "ImmuneCell") {
                agent->number = int(fe_para["frontend-api"]["immune_cell"]["number"].value("value", 10));
                agent->morphology_parameters.radius = fe_para["frontend-api"]["immune_cell"]["radius"].value("value", 10.6);
                agent->movement_parameters.mean = fe_para["frontend-api"]["immune_cell"]["speed"].value("value", 4);
                agent->movement_parameters.persistence_time = fe_para["frontend-api"]["immune_cell"]["persistence_time"].value("value", 1);
                bool diff_colors = fe_para["frontend-api"]["immune_cell"]["use_different_colors"].value("value", false);
                agent->morphology_parameters.color = diff_colors ? "random" : "greenTransp";
            }
        }
        sim_para.visualizer_to_overwrite = visp;
    }
}
