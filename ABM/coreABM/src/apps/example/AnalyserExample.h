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

#ifndef COREABM_ANALYSEREXAMPLE_H
#define COREABM_ANALYSEREXAMPLE_H

#include "InSituMeasurementsExample.h"
#include "core/analyser/Analyser.h"

class AnalyserExample : public Analyser{
  public:
    AnalyserExample(const std::string &config_path, const std::string &project_dir);
    std::shared_ptr<InSituMeasurements> generateMeasurement(const std::string &id) const override;

  protected:
    mutable std::vector<std::shared_ptr<InSituMeasurementsExample>> measurements_;
};

#endif // COREABM_ANALYSEREXAMPLE_H
