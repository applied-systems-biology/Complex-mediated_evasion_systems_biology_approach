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

#ifndef COREABM_ANALYSERCME_H
#define COREABM_ANALYSERCME_H

#include "core/analyser/Analyser.h"

class InSituMeasurementsCME;

class AnalyserCME : public Analyser {
  public:
    AnalyserCME(const std::string &config_path, const std::string &project_dir);
    ~AnalyserCME() override= default;
    std::shared_ptr<InSituMeasurements> generateMeasurement(const std::string &id) const override;
    void outputAllMeasurements() const override;

  protected:
    mutable std::vector<std::shared_ptr<InSituMeasurementsCME>> measurements_;
};

#endif // COREABM_ANALYSERCME_H
