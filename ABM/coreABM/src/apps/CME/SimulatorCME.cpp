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
#include "SimulatorCME.h"
#include <omp.h>
#include "apps/CME/cells/Complex.h"

#include "AnalyserCME.h"
#include "CuboidSiteCME.h"

void SimulatorCME::executeRuns(int runs, int seed, const std::string& output_dir, const std::string& input_dir, int sim,
                                   const std::string& parameter_string) const{
    SYSTEM_STDOUT("SimulatorCME engaged. Hang tight as we dive into the molecular world and simulate complex-mediated evasion of antimicrobial peptides!");

    Simulator::executeRuns(runs, seed, output_dir, input_dir, sim, parameter_string);
}

std::unique_ptr<Site>
SimulatorCME::createSites(int run, Randomizer* random_generator,
                              const Analyser* analyser) const {
    return std::make_unique<CuboidSiteCME>(random_generator,
                                        analyser->generateMeasurement(std::to_string(run)),
                                        config_path_, cmd_input_args_, output_dir_);

}

std::unique_ptr<const Analyser> SimulatorCME::createAnalyser(std::string config_path_, std::string project_dir) const {
    return std::make_unique<const AnalyserCME>(config_path_, project_dir);
}