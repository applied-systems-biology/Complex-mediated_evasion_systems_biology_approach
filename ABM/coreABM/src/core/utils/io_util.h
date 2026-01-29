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

#ifndef CORE_UTILS_IO_UTIL_H
#define CORE_UTILS_IO_UTIL_H

#include <iostream>
#include <string>
#include <unordered_set>
#include <map>
#include <unordered_map>
#include <array>
#include <type_traits>

#include <boost/filesystem.hpp>

#include "core/basic/Coordinate3D.h"

namespace abm::util {

    struct ConfigParameters {
        int runs{};
        int number_of_threads{};
        int system_seed{};
        std::string simulator{};
        std::string config_path{};
        std::string output_dir{};
        std::string input_dir{};
        std::unordered_map<std::string, std::vector<std::string>> screening_parameters{};
        int screen_start_idx{};
        int screen_stop_idx{};
    };

    struct VisualizerParameters {
        int run_id{};
        int total_runs{};
        int output_interval{};
        bool pov_active{};
        bool image_noise{};
        bool include_time{};
        bool output_video{};
        double white_noise{};
        double camera_angle{};
        std::string px_width{};
        std::string px_height{};
        std::string path_style_transfer{};
        Coordinate3D camera_position{};
        Coordinate3D camera_look_at{};
        std::vector<Coordinate3D> light_sources;
    };

    struct SimulationParameters {
        struct MovementParameters {
            double diffusion_coefficient{};
            double persistence_time{};
            double mean{};
            double stddev{};
            std::string type{};
        };
        struct MorphologyParameters {
            double radius{};
            double stddev{};
            std::string type{};
            std::string color{};
        };
        struct CellStateParameters {
            std::string name{};
            std::vector<std::pair<std::string, std::string>> next_states;
        };

        struct BinomialDistribution {
            bool activated{};
            std::uint64_t n{};
            double p{};
        };

        struct AgentParameters {
            int initial_distribution{};
            double input_lambda{};
            float number{};
            std::string input_distribution_path{};
            std::string type{};
            BinomialDistribution binomial_distribution{};
            MorphologyParameters morphology_parameters{};
            MovementParameters movement_parameters{};
            MovementParameters passive_movement_parameters{};
            std::vector<CellStateParameters> states;
        };

        struct ImmuneCellParameters : public AgentParameters {
        };

        struct FungalParameters : public AgentParameters {
        };

        struct AgentManagerParameters {
            std::string site_identifier{};
            std::vector<std::shared_ptr<AgentParameters>> agents;
        };

        struct NHLParameters {
            std::string type{};
            int interaction_check_interval{};
        };
        struct SiteParameters {
            std::string type{};
            std::string identifier{};
            std::string surface{};
            bool passive_movement{};
            std::string boundary_condition{};
            NHLParameters nhl_parameters{};
            AgentManagerParameters agent_manager_parameters{};
        };
        struct CuboidSiteParameters : public SiteParameters {
            Coordinate3D upper_bound{};
            Coordinate3D lower_bound{};
            std::tuple<int, int, int> molecules_grid_size{};
            //std::array<double> grid_size{};
        };
        struct InteractionStateParameters {
            bool adhere{};
            double must_overhead{};
            std::string name{};
            std::string interaction_type{};
            std::vector<std::pair<std::string, std::string>> next_states;
        };

        struct InteractionParameters {
            //name and vector for all states the condition is restricted to, if the vector is empty all states of all cell are ok for the condition
            std::string name{};
            std::string type{};
            std::vector<std::pair<std::string, std::vector<std::string>>> cell_conditions;
            std::vector<InteractionStateParameters> states{};
        };

        int dimensions{};
        bool use_interactions{};
        double max_time{};
        double time_stepping{};
        std::vector<std::string> stopping_criteria{};
        std::string topic{};
        std::vector<std::unique_ptr<InteractionParameters>> interaction_parameters;
        std::unique_ptr<SiteParameters> site_parameters;
        std::unordered_map<std::string, std::string> cmd_input_args{};
        VisualizerParameters visualizer_to_overwrite{};
    };

struct AnalyserParameters {
        std::unordered_set<std::string> active_measurements{};
        std::vector<std::string> cell_state_count{};
    };

    struct InputParameters {
        struct DefaultRateParameters {
            double rate{};
            std::string type{};
            std::string key{};
            // for conditional rate
            std::string condition{};
            virtual ~DefaultRateParameters() = default;
        };
        struct DistributionParameters {
            std::string type{};
            std::string key{};
            std::string source_file{};
            int number_of_distributions{};
            double lambda{};
        };
        std::vector<std::unique_ptr<DefaultRateParameters>> rates{};
        std::vector<DistributionParameters> distributions{};
    };

    InputParameters getInputParameters(const std::string &input_config);
    VisualizerParameters getViualizerParameters(const std::string &visualizer_config);
    AnalyserParameters getAnalyserParameters(const std::string &analyser_config);
    ConfigParameters getMainConfigParameters(const std::string &config_path);
    SimulationParameters getSimulationParameters(const std::string &simulator_config);

    void executeShellCommand(const std::string &command, bool suppress_output = true);

}
#endif //CORE_UTILS_IO_UTIL_H
