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
// Created by ybachelot on 24.04.23.
//

#include "InSituMeasurementsExample.h"
#include "core/simulation/Site.h"


InSituMeasurementsExample::InSituMeasurementsExample(std::unordered_set<std::string> active_measurements, const std::string& id){

    active_measurements_ = active_measurements;
    for (const auto& active_ : active_measurements_) {
        // Remove everything after "%" in the active measurement string
        auto active = active_;
        int char_pos = active_.find("%");
        if (char_pos > 0) {
            active = active_.substr(0, char_pos);
        }

        // Here define your own measurements
        if ("example_measurement" == active) {
            pair_measurements_["example"] = std::make_unique<PairMeasurement>(id, "example");
        }
    }
}
void InSituMeasurementsExample::observeMeasurements(const SimulationTime& time) {
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

        //Define here all your measurements
        // Measurements ever X timestep
        if ("example_measurement" == active && do_measurement) {
            for (const auto agent : site_->getAgentManager()->getAllAgents()) {
                pair_measurements_["example"]->addValuePairs(current_time,
                                                                      agent->getTypeName());
            }
        }
    }
}