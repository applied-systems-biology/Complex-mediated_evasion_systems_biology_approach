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

#ifndef CORE_ANALYSER_ANALYSER_H
#define CORE_ANALYSER_ANALYSER_H

#include <boost/filesystem.hpp>
#include <memory>

#include "core/utils/io_util.h"

class InSituMeasurements;

class Analyser {

public:
  // Class for wrapping analyzer functionality to take measurements during a simulation run
    Analyser() = default;
    virtual ~Analyser() = default;
    Analyser(const std::string &config_path, const std::string &project_dir);
    Analyser(const Analyser &) = delete;
    Analyser &operator=(const Analyser &) = delete;
    Analyser(Analyser &&) = delete;
    Analyser &operator=(Analyser &&) = delete;

    virtual std::shared_ptr<InSituMeasurements> generateMeasurement(const std::string &id) const;
    virtual void outputAllMeasurements() const;
protected:
    mutable std::vector<std::shared_ptr<InSituMeasurements>> measurements_;
    std::string measurement_path_;
    abm::util::AnalyserParameters parameters_;
};


#endif /* CORE_ANALYSER_ANALYSER_H */
