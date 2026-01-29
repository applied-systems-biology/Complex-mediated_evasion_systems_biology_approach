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

#ifndef COREABM_INSITUMEASUREMENTSCME_H
#define COREABM_INSITUMEASUREMENTSCME_H


#include "core/analyser/InSituMeasurements.h"


class InSituMeasurementsCME : public InSituMeasurements {
  public:
    InSituMeasurementsCME(std::unordered_set<std::string> active_measurements, const std::string &);
    ~InSituMeasurementsCME() override = default;
    void observeMeasurements(const SimulationTime &time);
    void writeToFiles(const std::string& output_dir) const override;


};

#endif // COREABM_INSITUMEASUREMENTSCME_H
