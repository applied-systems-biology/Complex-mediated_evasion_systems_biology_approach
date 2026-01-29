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
// Created by ybachelot on 21.04.23.
//

#include "AnalyserCME.h"
#include "apps/CME/InSituMeasurementsCME.h"
#include "core/analyser/Analyser.h"

AnalyserCME::AnalyserCME(const std::string& config_path, const std::string& project_dir)
    : Analyser( std::move(config_path), std::move(project_dir)){}


std::shared_ptr<InSituMeasurements> AnalyserCME::generateMeasurement(const std::string& id) const {
    std::shared_ptr<InSituMeasurementsCME> new_measurement = nullptr;
#pragma omp critical
    {
        new_measurement = measurements_.emplace_back(
            std::make_shared<InSituMeasurementsCME>(parameters_.active_measurements, id));
    }
    return new_measurement;
}

void AnalyserCME::outputAllMeasurements() const {

    for (const auto &measurement: measurements_) {
        measurement->writeToFiles(measurement_path_);
    }
}