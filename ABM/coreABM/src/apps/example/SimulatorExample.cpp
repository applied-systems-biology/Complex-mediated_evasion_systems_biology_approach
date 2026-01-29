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

#include "core/simulation/Site.h"
#include "core/utils/macros.h"
#include "SimulatorExample.h"

#include "AnalyserExample.h"
#include "CuboidSiteExample.h"
#include <unordered_map>

void SimulatorExample::executeRuns(int runs, int seed, const std::string& output_dir, const std::string& input_dir, int sim,
                                   const std::string& parameter_string) const {
    SYSTEM_STDOUT("SimulatorExample: This is going to be amazing...");
    Simulator::executeRuns(runs, seed, output_dir, input_dir, sim, parameter_string);
}
std::unique_ptr<Site>
SimulatorExample::createSites(int run, Randomizer* random_generator,
                              const Analyser* analyser) const {
    SYSTEM_STDOUT("SimulatorExample: This will even be more amazing...");

    return std::make_unique<CuboidSiteExample>(random_generator,
                                        analyser->generateMeasurement(std::to_string(run)),
                                        config_path_, cmd_input_args_, output_dir_);

}

std::unique_ptr<const Analyser> SimulatorExample::createAnalyser(std::string config_path_, std::string project_dir) const {
    return std::make_unique<const AnalyserExample>(config_path_, project_dir);
}
