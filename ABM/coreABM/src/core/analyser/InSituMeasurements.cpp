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

#include <cmath>

#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/Site.h"
#include "core/simulation/neighbourhood/Collision.h"
#include "core/simulation/Interaction.h"
#include "core/utils/macros.h"

InSituMeasurements::InSituMeasurements(std::unordered_set<std::string> active_measurements, const std::string &id)
        : active_measurements_(std::move(active_measurements)) {
    for (const auto &active_:active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            active = active_.substr(0, char_pos);
        }

        if ("agent-statistics" == active) {
            time_last_measurement_["agent-statistics"] = 0.0;
            pair_measurements_["agent-statistics"] = std::make_unique<PairMeasurement>(id, "time", "agent", "agentid",
                                                                                       "state", "radius", "x", "y", "z",
                                                                                       "cellpart", "cellpart_id");
        } else if ("environment" == active) {
            pair_measurements_["environment"] = std::make_unique<PairMeasurement>(id, "x", "y", "z", "radius_or_length",
                                                                                  "type", "additional");
        }
    }
}

void InSituMeasurements::observeMeasurements(const SimulationTime &time) {
    using std::make_pair;
    const auto &current_time = time.getCurrentTime();

    for (const auto &active_:active_measurements_) {

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
            for (const auto agent: site_->getAgentManager()->getAllAgents()) {
                std::string cellpart = "Mothercell";
                for (const auto cellparts: agent->getMorphology()->getAllSpheresOfThis()) {
                    if (cellparts->getCreationTime() >= time_last_measurement_["agent-statistics"] || cellpart != "Hyphae") {
                        pair_measurements_["agent-statistics"]->addValuePairs(current_time,
                                                                              agent->getTypeName(), agent->getId(),
                                                                              agent->getCurrentCellState()->getStateName(),
                                                                              cellparts->getRadius(),
                                                                              cellparts->getPosition().x,
                                                                              cellparts->getPosition().y,
                                                                              cellparts->getPosition().z,
                                                                              cellpart, cellparts->getDescription());
                        // First cellpart is always "Mothercell" - everything else is Hyphae
                        cellpart = "Hyphae";
                    }
                }
            }
            time_last_measurement_["agent-statistics"] = current_time;
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
        }
    }
}

void InSituMeasurements::writeToFiles(const std::string &output_dir) const {
    for (const auto&[name, measurement]: histogram_measurements_) {
        const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
        if (!boost::filesystem::exists(file_name)) {
            std::ofstream file{file_name};

            const auto &keys = measurement->getKeys();
            file << "id" << HistogramMeasurement::delimeter;
            for (const auto &key: keys) {
                file << key << HistogramMeasurement::delimeter;
            }
            file << '\n';
            file << *measurement;
            file.close();
            measurement->clearData();
        } else {
            std::ofstream file{file_name, std::ios_base::app};
            file << *measurement;
            file.close();
            measurement->clearData();
        }
    }
    for (const auto&[name, measurement]: pair_measurements_) {
        const auto file_name = static_cast<boost::filesystem::path>(output_dir).append(name + ".csv").string();
        if (!boost::filesystem::exists(file_name)) {
            std::ofstream file{file_name};
            const auto &keys = measurement->getKeys();
            file << "id" << PairMeasurement::delimeter;
            for (const auto &key: keys) {
                file << key << PairMeasurement::delimeter;
            }
            file << '\n';
            file << *measurement;
            file.close();
            measurement->clearData();
        } else {
            std::ofstream file{file_name, std::ios_base::app};
            file << *measurement;
            file.close();
            measurement->clearData();
        }
    }
}