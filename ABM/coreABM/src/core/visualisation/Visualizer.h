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

#ifndef CORE_VISUALISATION_VISUALIZER_H
#define CORE_VISUALISATION_VISUALIZER_H

#include <string>
#include "core/utils/io_util.h"
#include "core/utils/time_util.h"


class Site;

class Visualizer {
public:
    /// Central class for visualizing a simulation
    Visualizer() = default;
    Visualizer(const std::string &config_path, const std::string &project_dir, int total_runs = 0);

    /*!
     * Visualize current configuration per timestep
     * @param site Site object of environment (e.g. AlveoleSite)
     * @param time SimulationTime object, i.e. contains current time and timestep
     * @param run Integer that contains current run
     * @param simulation_end Boolean that denotes if simulation end is reached
     */
    virtual void visualizeCurrentConfiguration(Site &site, const SimulationTime &time, int run,
                                       bool simulation_end = false) const;

    virtual void overwriteParameters(abm::util::VisualizerParameters new_parameters) const {
        if (new_parameters.output_interval != 0) parameters_ = new_parameters;};

protected:
    void concludeRun(Site &site, int run) const;
    std::vector<std::string> visualization_path_{};
    mutable abm::util::VisualizerParameters parameters_{};
};

#endif /* CORE_VISUALISATION_VISUALIZER_H */
