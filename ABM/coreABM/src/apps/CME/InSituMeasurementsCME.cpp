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
// Created by ybachelot on 19.04.23.
//

#include "InSituMeasurementsCME.h"
#include <cmath>

#include "apps/CME/cells/Complex.h"
#include "apps/CME/cells/Pathogenic_cell.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/Interaction.h"
#include "core/simulation/Site.h"
#include "core/simulation/neighbourhood/Collision.h"
#include "core/utils/macros.h"

InSituMeasurementsCME::InSituMeasurementsCME(std::unordered_set<std::string> active_measurements, const std::string& id){

    active_measurements_ = active_measurements;
    for (const auto& active_ : active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            active = active_.substr(0, char_pos);
        }

        if ("agent-statistics" == active) {
            pair_measurements_["agent-statistics"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                       "state", "radius", "x", "y", "z",
                                                                                       "dist_to_cell");
        } else if ("environment" == active) {
            pair_measurements_["environment"] = std::make_unique<PairMeasurement>(id, "x", "y", "z", "radius_or_length",
                                                                                  "type", "additional");
        } else if ("lifetime-distribution" == active) {
            histogram_measurements_["lifetime"] = std::make_unique<HistogramMeasurement>(id, "timestep", "lifetime_complex");
        } else if ("bounds-distrib" == active) {
            histogram_measurements_["bounds_distrib"] = std::make_unique<HistogramMeasurement>(id, "agent", "nb_bounds");
        } else if ("bounds-uptake" == active) {
            histogram_measurements_["bounds_uptake"] = std::make_unique<HistogramMeasurement>(id, "nb_bounds");
        } else if ("total-uptake" == active) {
            pair_measurements_["uptake"] = std::make_unique<PairMeasurement>(id, "time", "uptake");
        } else if ("population-statistics" == active){
            pair_measurements_["pop"] = std::make_unique<PairMeasurement>(id, "time", "agent", "nb");
        } else if ("agent-statistics-multiple_receptors" == active) {
            pair_measurements_["agent-statistics-multiple_receptors"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                       "state", "radius", "x", "y", "z",
                                                                                       "dist_to_cell", "nb_AMP");
        } else if ("agent-statistics-multiple_receptors_drug" == active) {
            pair_measurements_["agent-statistics-multiple_receptors_drug"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                                          "state", "radius", "x", "y", "z",
                                                                                                          "dist_to_cell", "nb_AMP", "nb_Drug");
        } else if ("proba_survival" == active){
            pair_measurements_["proba_survival"] = std::make_unique<PairMeasurement>(id, "time", "proba_survival");
        }
        else if ("uptake_rate" == active) {
            pair_measurements_["uptake_rate"] = std::make_unique<PairMeasurement>(id, "time", "uptake_rate");
        }
        else if ("complex-content" == active) {
            pair_measurements_["complex-content"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                       "state", "radius", "nb_AMP", "nb_drug", "free_receptor");
        }
        else if ("complex_formation" == active){
            histogram_measurements_["complex_formation"] = std::make_unique<HistogramMeasurement>(id, "time", "molecule");
        }
    }
}

void InSituMeasurementsCME::observeMeasurements(const SimulationTime &time){
    using std::make_pair;
    const auto& current_time = time.getCurrentTime();

    for (const auto& active_ : active_measurements_) {

        // Code snippet to let your measurement only execute ever X-th simulation step
        // Name your measurement as "measurment%X" -> example: "agent-statistics%25" ... every 25th step
        // Default if no value is given: Every 10th step
        int every_x_step = 10;
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            every_x_step = stoi(active_.substr(char_pos + 1, active_.size()));
            active = active_.substr(0, char_pos);
        }
        bool do_measurement = ((time.getCurrentTimeStep() % every_x_step) == 0) || site_->checkForStopping() || time.lastStepBeforEnd();

        // Measurements ever X timestep
        if ("agent-statistics" == active && do_measurement) {
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                pair_measurements_["agent-statistics"]->addValuePairs(current_time,
                                                                      agent->getTypeName(), agent->getId(),
                                                                      agent->getCurrentCellState()->getStateName(),
                                                                      agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                      agent->getPosition().x,
                                                                      agent->getPosition().y,
                                                                      agent->getPosition().z,
                                                                      agent->getPosition().calculateEuclidianDistance(Coordinate3D{0, 0, 0}));
            }
        }
        if ("agent-statistics-multiple_receptors" == active && do_measurement) {
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                if (agent->getTypeName() != "Complex"){
                    pair_measurements_["agent-statistics-multiple_receptors"]->addValuePairs(current_time,
                                                                                             agent->getTypeName(), agent->getId(),
                                                                                             agent->getCurrentCellState()->getStateName(),
                                                                                             agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                                             agent->getPosition().x,
                                                                                             agent->getPosition().y,
                                                                                             agent->getPosition().z,
                                                                                             agent->getPosition().calculateEuclidianDistance(Coordinate3D{0, 0, 0}),
                                                                                             0);
                }
                else {
                    auto nb = std::static_pointer_cast<Complex>(agent)->get_complex_content().size();
                    pair_measurements_["agent-statistics-multiple_receptors"]->addValuePairs(current_time,
                                                                                             agent->getTypeName(), agent->getId(),
                                                                                             agent->getCurrentCellState()->getStateName(),
                                                                                             agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                                             agent->getPosition().x,
                                                                                             agent->getPosition().y,
                                                                                             agent->getPosition().z,
                                                                                             agent->getPosition().calculateEuclidianDistance(Coordinate3D{0, 0, 0}),
                                                                                             nb);
                }
            }
        }
        if ("agent-statistics-multiple_receptors_drug" == active && do_measurement) {
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                if (agent->getTypeName() != "Complex"){
                    pair_measurements_["agent-statistics-multiple_receptors_drug"]->addValuePairs(current_time,
                                                                                             agent->getTypeName(), agent->getId(),
                                                                                             agent->getCurrentCellState()->getStateName(),
                                                                                             agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                                             agent->getPosition().x,
                                                                                             agent->getPosition().y,
                                                                                             agent->getPosition().z,
                                                                                             agent->getPosition().calculateEuclidianDistance(Coordinate3D{0, 0, 0}),
                                                                                             0,
                                                                                             0);
                }
                else {
                    auto complex_content = std::static_pointer_cast<Complex>(agent)->get_complex_content();
                    auto nb_AMP = std::count(complex_content.begin(), complex_content.end(), "AMP");;
                    auto nb_drug = std::count(complex_content.begin(), complex_content.end(), "Drug");

                    pair_measurements_["agent-statistics-multiple_receptors_drug"]->addValuePairs(current_time,
                                                                                             agent->getTypeName(), agent->getId(),
                                                                                             agent->getCurrentCellState()->getStateName(),
                                                                                             agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                                             agent->getPosition().x,
                                                                                             agent->getPosition().y,
                                                                                             agent->getPosition().z,
                                                                                             agent->getPosition().calculateEuclidianDistance(Coordinate3D{0, 0, 0}),
                                                                                             nb_AMP,
                                                                                             nb_drug);
                }
            }
        }
        if ("complex-content" == active && do_measurement) {
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                if (agent->getTypeName() == "Complex"){
                    auto complex_content = std::static_pointer_cast<Complex>(agent)->get_complex_content();
                    auto nb_AMP = std::count(complex_content.begin(), complex_content.end(), "AMP");
                    auto nb_drug = std::count(complex_content.begin(), complex_content.end(), "Drug");
                    auto free_recep = std::static_pointer_cast<Complex>(agent)->get_free_receptors();
                    pair_measurements_["complex-content"]->addValuePairs(current_time,agent->getTypeName(), agent->getId(),
                                                                                      agent->getCurrentCellState()->getStateName(),
                                                                                      agent->getMorphology()->getBasicSphereOfThis()->getRadius(),
                                                                                      nb_AMP,
                                                                                      nb_drug,
                                                                                      free_recep);
                }
            }
        }
        if ("total-uptake" == active && do_measurement) {
            const auto bounds = site_->getSystemBoundaries();
            const auto volume = (bounds[1].x - bounds[0].x) *(bounds[1].y - bounds[0].y)*(bounds[1].z - bounds[0].z);
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                if (agent->getTypeName() == "Pathogenic_cell") {
                    const auto ag = std::dynamic_pointer_cast<Pathogenic_cell>(agent);
                    pair_measurements_["uptake"]->addValuePairs(current_time, ag->get_total_uptake()/volume);
                }
            }
        }
        if ("uptake_rate" == active && do_measurement) {
            pair_measurements_["uptake_rate"]->addValuePairs(current_time, site_->getRateFactory()->getRate("uptake")->getRateValue());
        }
        if ("population-statistics" == active && do_measurement){
            const auto bounds = site_->getSystemBoundaries();
            const auto volume = (bounds[1].x - bounds[0].x) *(bounds[1].y - bounds[0].y)*(bounds[1].z - bounds[0].z);
            for (const auto& name : site_->getAgentManager()->getAllAgentTypes()) {
                const auto nb = site_->getAgentManager()->getAgentQuantity(name);
                pair_measurements_["pop"]->addValuePairs(current_time, name, nb/volume);
            }
        }
        if ("proba_survival" == active && do_measurement){
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                if (agent->getTypeName() == "Pathogenic_cell") {
                    const auto ag = std::dynamic_pointer_cast<Pathogenic_cell>(agent);
                    pair_measurements_["proba_survival"]->addValuePairs(current_time, ag->get_proba_survival());
                }
            }
        }

        // End of the simulation
        if (time.getMaxTime() - time.getCurrentTime() <= time.getCurrentDeltaT() or site_->checkForStopping()) {
            if ("environment" == active) {
                std::vector<Coordinate3D> sysBound = site_->getSystemBoundaries();
                Coordinate3D pointA = sysBound[0];
                Coordinate3D pointB = sysBound[1];
                pair_measurements_["environment"]->addValuePairs(pointA.x, pointA.y, pointA.z, "0.0",
                                                                 "minBounds", 0.0);
                pair_measurements_["environment"]->addValuePairs(pointB.x, pointB.y, pointB.z, "0.0",
                                                                 "maxBounds", 0.0);
            }

            else if ("bounds-distrib" == active){
                /*
                for (auto const &ag: site_->getAgentManager()->getAllAgents()) {
                    if (ag->getTypeName() == "AMP" || ag->getTypeName() == "Defensive"){
                        const auto cell = std::dynamic_pointer_cast<Cell>(ag);
                        auto agt = 0;
                        if (cell->getTypeName() == "AMP"){
                            agt = 1;
                        }
                        else{
                            agt = 2;
                        }
                        addValues<HistogramMeasurement>("bounds_distrib", std::make_pair("agent", agt), std::make_pair("nb_bounds", cell->get_nb_bound()));
                    }
                    else if (ag->getTypeName() == "Complex"){
                        const auto cell = std::dynamic_pointer_cast<Complex>(ag);
                        auto bounds = cell->get_complex_content().size();
                        for (auto bound: std::get<0>(bounds)){
                            addValues<HistogramMeasurement>("bounds_distrib", std::make_pair("agent", 1), std::make_pair("nb_bounds", bound));
                        }
                        addValues<HistogramMeasurement>("bounds_distrib", std::make_pair("agent", 2), std::make_pair("nb_bounds", std::get<1>(bounds)));
                    }
                }
                 */
            }
            // SAVE LIFETIME STATS
            for (auto const ag: site_->getAgentManager()->getAllAgents()) {
                if (ag->getTypeName() == "Complex") {
                    auto lifetime = std::dynamic_pointer_cast<Complex>(ag)->get_lifetime();
                    addValues<HistogramMeasurement>("lifetime", std::make_pair("timestep", time.getCurrentDeltaT()), std::make_pair("lifetime_complex", lifetime));
                }
            }
        }
    }
}

void InSituMeasurementsCME::writeToFiles(const std::string& output_dir) const {
    InSituMeasurements::writeToFiles(output_dir);

    // Generate HTML files containing plots
    std::string command = "/home/ybachelot/miniconda3/bin/conda"; // conda folder from json?
    command.append(" run -n abmenv python ");
    command.append("src/apps/CME/python_scripts/output_graphs.py ");
    command.append(output_dir);

    system(command.c_str());
/*
    // Generate xyz files for visualization in ovito
    command = "/home/ybachelot/miniconda3/bin/conda";
    command.append(" run -n abmenv python ");
    command.append("src/apps/CME/python_scripts/ovito_preprocessing.py ");
    command.append(output_dir);
    system(command.c_str());
*/
}