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

#include <boost/filesystem.hpp>
#include <iostream>
#include <map>

#include "apps/example/SimulatorExample.h"
#include "apps/CME/SimulatorCME.h"
#include "core/simulation/Simulator.h"
#include "core/utils/io_util.h"
#include "core/utils/macros.h"
#include "core/utils/misc_util.h"
#include <omp.h>


std::unique_ptr<Simulator> createSimulator(const std::string simulator) {
    // Different Simulators are defined for different subprojects. It takes the value for "simulator" in the main config.json.
    // Please write a short description when you add a new Simulator.
    // -> "Simulator" is the simulator for the coreABM and only uses functionalities in the "core/" folder
    // -> "SimulatorExample" is an example that can be used to create a new simulator and is located in the "simulator/" folder
    if (simulator == "SimulatorExample") {
        return std::make_unique<SimulatorExample>();
    } else if (simulator == "SimulatorCME"){
        return std::make_unique<SimulatorCME>();
    } else {
        return std::make_unique<Simulator>();
    }
}

int main(int argc, char** argv) {

    // Default location of config file (if no parameter is specified)
    boost::filesystem::path config_xml("../../config.json");
    if (argc > 1 && !(std::istringstream(argv[1]) >> config_xml)) {
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 1;
    }
    if (!(boost::filesystem::exists(config_xml))) {
        ERROR_STDERR("Configuration File in " << config_xml.string() << " does not exist!");
        ERROR_STDERR("usage: " << argv[0] << " <config.json>");
        return 2;
    }

    // Change root directory for simulation to configuration location
    auto chd = chdir(&config_xml.parent_path().c_str()[0]);

    // Read cmd inputs and screening parameters
    auto parameters = abm::util::getMainConfigParameters(config_xml.filename().string());
    std::unordered_map<std::string, std::string> input_args = abm::util::handleCmdInputs(argc, argv);

    // Initialize parallelization
#if defined(_OPENMP)
    int max_threads_available = omp_get_max_threads();
    // Set the number of threads to the minimum of parameters.number_of_threads and max_threads_available
    int nb_threads = std::min(parameters.number_of_threads, max_threads_available);
    omp_set_num_threads(nb_threads);
    omp_set_nested(1); // Enable nested parallelism
    DEBUG_STDOUT("OpenMP activated with " << parameters.number_of_threads << " Thread(s).");
#endif

    // Start simulation runs
    if (parameters.screening_parameters.empty()) {
//        const auto simulator = std::make_unique<const SimulatorExample>(parameters.config_path, input_args);
        auto simulator = createSimulator(parameters.simulator);
        simulator->setConfigPath(parameters.config_path);
        simulator->setCmdInputArgs(input_args);
        simulator->setOutputPath(parameters.output_dir);
        simulator.get()->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir, parameters.input_dir);
        simulator.reset();
    } else {
        // Screening over all parameter combinations specified as sets in the configuration file <config.json>
        // For screening, the cartesian product of all the single parameters sets is generated
        // If you want to resume a screening from a certain index, you can change screen_start_idx in main config
        auto result = abm::util::calculateCartesianProd(parameters.screening_parameters);
        const auto& parameter_names = result.first;
        const auto& value_combinations = result.second;
        if (parameters.screen_stop_idx == 0){
            parameters.screen_stop_idx = value_combinations.size();
        }

        std::vector<std::unordered_map<std::string, std::string>> local_input_args_list(parameters.number_of_threads); // Vector to store intermediate results

        int nb_threads_to_use = 1;
        int rep = parameters.runs;
        if (nb_threads > rep){
            nb_threads_to_use = std::floor(nb_threads / rep);
        }
#pragma omp parallel for schedule(dynamic) num_threads(nb_threads_to_use)
        for (int sim = parameters.screen_start_idx-1; sim < parameters.screen_stop_idx; ++sim) {
            int thread_id = omp_get_thread_num();
            std::unordered_map<std::string, std::string> local_input_args;
            std::stringstream sim_para{};
            for (int i = 0; i < parameter_names.size(); ++i) {
                sim_para << parameter_names[i] << value_combinations[sim][i] << "_";
                local_input_args[parameter_names[i]] = value_combinations[sim][i];
            }
            local_input_args_list[thread_id] = std::move(local_input_args);
            SYSTEM_STDOUT("Start " << parameters.runs << " runs of simulation " << sim + 1 << "/" << parameters.screen_stop_idx);
//            auto simulator = std::make_unique<const SimulatorExample>(parameters.config_path, input_args);
            auto simulator = createSimulator(parameters.simulator);
            simulator->setConfigPath(parameters.config_path);
            simulator->setCmdInputArgs(local_input_args_list[thread_id]);
            simulator->setOutputPath(parameters.output_dir);
            simulator.get()->executeRuns(parameters.runs, parameters.system_seed, parameters.output_dir, parameters.input_dir, sim, sim_para.str());
            simulator.reset();
        }
    }
    return 0;
}
