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

#ifndef CORE_SIMULATION_SIMULATOR_H
#define CORE_SIMULATION_SIMULATOR_H

#include "core/utils/io_util.h"
#include "core/utils/time_util.h"
#include "core/visualisation/Visualizer.h"

class Site;
class Analyser;
class Randomizer;

namespace abm::test { std::string test_simulation(const std::string &config); }
class Simulator {
public:
    /// Class for starting simulations
    Simulator() {};
    Simulator(std::string config_path, std::unordered_map<std::string, std::string> cmd_input_args);
    ~Simulator() = default;
    Simulator(const Simulator &) = delete;
    auto operator=(const Simulator &) -> Simulator & = delete;
    Simulator(Simulator &&) = delete;
    auto operator=(Simulator &&) -> Simulator & = delete;

    /**
     * Initializes all objects needed for a simulation and contains the main for-loop over all timesteps
     * @param runs Integer that contains the number of runs
     * @param seed Integer that contains the seed value
     * @param output_dir String that contains the output directory
     * @param input_dir String that contains the input directory
     * If parameter screening is activated:
     * @param sim Integer that contains the enumerated parameter configuration that is currently used
     * @param parameter_string String that contains the parameter configuration that is currently used
     */
    virtual void executeRuns(int runs, int seed, const std::string &output_dir, const std::string &input_dir, int sim = 0,
                     const std::string &parameter_string = "") const;
    /**
     *  Creates visualizer for simulation
     * @param config_path_ path to configuration
     * @param project_dir name of project folder where results are stored
     * @param runs number of runs
     * @return an instance of class Visualizer
     */
    virtual std::unique_ptr<const Visualizer> createVisualizer(std::string config_path_, std::string project_dir, int runs) const;

    /**
     * Creates analyer for simulation to write outputs in "measurement" folder
     * @param config_path_ path to configuration
     * @param project_dir name of project folder where results are stored
     * @return an instance of class Analyser
     */
    virtual std::unique_ptr<const Analyser> createAnalyser(std::string config_path_, std::string project_dir) const;

    /*!
     * Creates the environment for the simulation
     * @param run Integer that contains the current run
     * @param random_generator Randomizer object
     * @param analyser Analyser object
     * @param input_dir String that contains input directory
     * @return Object of created Site (e.g. AlveoleSite)
     */
    virtual std::unique_ptr<Site> createSites(int run,
                                      Randomizer *random_generator,
                                      const Analyser *analyser) const;

    /// Functions to write output.json which is necessary for web-frontend gui
    void setConfigPath(std::string config_path) {config_path_ = config_path;}
    void setCmdInputArgs(std::unordered_map<std::string, std::string> cmd_input_args) {cmd_input_args_ = cmd_input_args;}
    void setOutputPath(std::string output_path) {output_dir_ = output_path;}
    void initFrontendAPIOutput(std::vector<std::pair<std::string, std::vector<std::string>>>& output, const std::string &project_name,
                               const int runs, const Analyser* analyser) const;
    void writeFrontendAPIOutput(std::vector<std::pair<std::string, std::vector<std::string>>>& output, Site &site, SimulationTime time,
                                const std::string &project_name, const int run) const;
    /// Used for integration tests
    friend std::string abm::test::test_simulation(const std::string &config);

protected:
    std::string config_path_{};
    std::string output_dir_{};
    std::unordered_map<std::string, std::string> cmd_input_args_{};
};

#endif // CORE_SIMULATION_SIMULATOR_H
