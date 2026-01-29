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

#include <optional>
#include <random>

#include "core/utils/io_util.h"
#include "external/json.hpp"
#include "core/utils/macros.h"

namespace abm::util {

    using json = nlohmann::json;

    ConfigParameters getMainConfigParameters(const std::string &main_config) {
        ConfigParameters parameters{};
        std::ifstream json_file(main_config);
        json json_parameters;
        json_file >> json_parameters;
        parameters.runs = json_parameters["Agent-Based-Framework"].value("runs", 1);
        parameters.simulator = json_parameters["Agent-Based-Framework"].value("simulator", "Simulator");
        parameters.config_path = json_parameters["Agent-Based-Framework"]["available_simulators"][parameters.simulator].value("config_path", "");
        parameters.output_dir = json_parameters["Agent-Based-Framework"]["available_simulators"][parameters.simulator].value("output_dir", "/tmp/results");
        parameters.input_dir = json_parameters["Agent-Based-Framework"].value("input_dir", "input/");
        parameters.number_of_threads = json_parameters["Agent-Based-Framework"].value("number_of_threads", 1);
        std::mt19937 gen(std::random_device{}());
        std::uniform_int_distribution dis;
        parameters.system_seed = json_parameters["Agent-Based-Framework"].value("seed", dis(gen));

        for (const auto &item: json_parameters["Agent-Based-Framework"]["parameter_screening"].items()) {
            std::vector<std::string> values;
            for (float value: item.value()) {
                std::string valuestr = std::to_string(value);
                if (std::stoi(valuestr.substr(valuestr.find_first_of('.') + 1, -1)) == 0) {
                    values.push_back(std::to_string(static_cast<int>(value)));
                } else {
                    values.push_back(std::to_string(value));
                }
            }
            parameters.screening_parameters[item.key()] = values;
        }
        parameters.screen_start_idx = json_parameters["Agent-Based-Framework"].value("screen_start_idx", 1);
        parameters.screen_stop_idx = json_parameters["Agent-Based-Framework"].value("screen_stop_idx", 0);

        json_file.close();
        return parameters;
    }

    SimulationParameters getSimulationParameters(const std::string &simulator_config) {
        SimulationParameters parameters{};
        std::ifstream json_file(simulator_config);
        json json_parameters;
        json_file >> json_parameters;
        try {
            parameters.topic = json_parameters["Agent-Based-Framework"]["topic"];
            parameters.dimensions = json_parameters["Agent-Based-Framework"]["dimensions"];
            parameters.max_time = json_parameters["Agent-Based-Framework"]["max_time"];
            parameters.time_stepping = json_parameters["Agent-Based-Framework"]["timestepping"];
            for (const auto &stopping_crit: json_parameters["Agent-Based-Framework"]["stopping_criteria"]) {
                parameters.stopping_criteria.emplace_back(stopping_crit);
            }
        } catch (const std::bad_optional_access &e) {
            ERROR_STDERR(e.what());
            throw;
        }
        for (const auto &site: json_parameters["Agent-Based-Framework"]["Sites"]) {
            auto site_para = std::make_unique<SimulationParameters::SiteParameters>();
            const auto type = site["type"];
            if ("CuboidSite" == type) {
                auto cs_para = std::make_unique<SimulationParameters::CuboidSiteParameters>();
                cs_para->lower_bound = {site["CuboidSite"]["x_range"][0], site["CuboidSite"]["y_range"][0], site["CuboidSite"]["z_range"][0]};
                cs_para->upper_bound = {site["CuboidSite"]["x_range"][1], site["CuboidSite"]["y_range"][1], site["CuboidSite"]["z_range"][1]};
                site_para = std::move(cs_para);
            }
            site_para->type = type;
            site_para->identifier = site["identifier"];
            site_para->boundary_condition = site.value("boundary_condition", "ReflectingBoundaries");

            // load neighbourhood locator
            site_para->nhl_parameters = {site["NeighbourhoodLocator"]["type"],
                                         site["NeighbourhoodLocator"].value("interaction_check_interval", 1)};

            // load agent manager
            // we need to keep an ordering for reproducing simulations, changing it leads to agents get differently initialized due to random values
            for (const auto &agent_type:site["AgentManager"]["Types"]) {
                std::shared_ptr<SimulationParameters::AgentParameters> agent_parameters{};
                auto agent = site["AgentManager"]["Agents"].at(std::string(agent_type));

                if (agent_type == "ImmuneCell") {
                    auto ic_parameters = std::make_shared<SimulationParameters::ImmuneCellParameters>();
                    agent_parameters = std::move(ic_parameters);

                } else if (agent_type == "FungalCell") {
                    auto fc_parameters = std::make_shared<SimulationParameters::FungalParameters>();
                    agent_parameters = std::move(fc_parameters);
                } else {
                    agent_parameters = std::make_shared<SimulationParameters::AgentParameters>();
                }
                agent_parameters->initial_distribution = agent.value("initial_distribution", 0);
                agent_parameters->input_lambda = agent.value("input_lambda", 0.0);
                agent_parameters->morphology_parameters.color = agent["Morphology"]["color"];

                agent_parameters->morphology_parameters.type = "SphericalMorphology";
                agent_parameters->morphology_parameters.radius = agent["Morphology"]["SphericalMorphology"]["radius"];
                agent_parameters->morphology_parameters.stddev = agent["Morphology"]["SphericalMorphology"]["stddev"];
                if (agent.find("Movement") != agent.end()) {
                    auto movement_type = agent["Movement"].value("type", "default");
                    agent_parameters->movement_parameters.type = movement_type;
                    if ("RandomWalk" == movement_type) {
                        agent_parameters->movement_parameters.diffusion_coefficient = agent["Movement"].value(
                                "diffusion_coefficient", 0.0);
                        if (agent["Movement"].find("speed") != agent["Movement"].end()) {
                            agent_parameters->movement_parameters.mean = agent["Movement"]["speed"].value("mean", 0.0);
                            agent_parameters->movement_parameters.stddev = agent["Movement"]["speed"].value("stddev",
                                                                                                            0.0);
                        }
                    } else if ("PersistentRandomWalk" == movement_type) {
                        agent_parameters->movement_parameters.persistence_time = agent["Movement"].value(
                                "persistence_time", 0.0);
                        if (agent["Movement"].find("speed") != agent["Movement"].end()) {
                            agent_parameters->movement_parameters.mean = agent["Movement"]["speed"].value("mean", 0.0);
                            agent_parameters->movement_parameters.stddev = agent["Movement"]["speed"].value("stddev",
                                                                                                            0.0);
                        }
                    }
                }
                if (agent.find("Passive Movement") != agent.end()) {
                    auto passive_movement_type = agent["Passive Movement"]["type"];
                    agent_parameters->passive_movement_parameters.type = passive_movement_type;
                    if ("RandomWalk" == passive_movement_type) {
                        agent_parameters->passive_movement_parameters.diffusion_coefficient =
                                agent["Passive Movement"].value("diffusion_coefficient", 0.0);
                        if (agent["Passive Movement"].find("speed") != agent["Passive Movement"].end()) {
                            agent_parameters->passive_movement_parameters.mean = agent["Passive Movement"]["speed"].value(
                                    "mean", 0.0);
                            agent_parameters->passive_movement_parameters.stddev = agent["Passive Movement"]["speed"].value(
                                    "stddev", 0.0);
                        }
                    } else if ("PersistentRandomWalk" == passive_movement_type) {
                        agent_parameters->passive_movement_parameters.persistence_time = agent["Passive Movement"].value(
                                "persistence_time", 0.0);
                        if (agent["Passive Movement"].find("speed") != agent["Passive Movement"].end()) {
                            agent_parameters->passive_movement_parameters.mean = agent["Passive Movement"]["speed"].value(
                                    "mean", 0.0);
                            agent_parameters->passive_movement_parameters.stddev = agent["Passive Movement"]["speed"].value(
                                    "stddev", 0.0);
                        }
                    }
                }
                agent_parameters->type = agent_type;
                agent_parameters->number = agent.value("number", 0);
                agent_parameters->binomial_distribution.activated = agent.value("binomial", false);

                for (const auto &state:agent["Cell States"].items()) {
                    SimulationParameters::CellStateParameters cell_state_parameters{};
                    cell_state_parameters.name = state.key();
                    if (state.value().find("next_states") != state.value().end()) {
                        for (const auto &next_state:state.value()["next_states"].items()) {
                            cell_state_parameters.next_states.emplace_back(next_state.key(),
                                                                           next_state.value()["rate"]);
                        }
                    }
                    agent_parameters->states.emplace_back(cell_state_parameters);
                }
                site_para->agent_manager_parameters.agents.emplace_back(std::move(agent_parameters));
            }


            parameters.site_parameters = std::move(site_para);
        }

        parameters.use_interactions = json_parameters["Agent-Based-Framework"].value("use_interactions", true);
        if (auto interactions = json_parameters["Agent-Based-Framework"].find("Interactions");interactions
                                                                                              !=
                                                                                              json_parameters["Agent-Based-Framework"].end()) {
            for (const auto&[key, interaction]:interactions->items()) {
                auto interaction_parameters = std::make_unique<SimulationParameters::InteractionParameters>();
                interaction_parameters->type = interaction["type"];
                interaction_parameters->name = key;
                interaction_parameters->states = {};
                for (const auto &state:interaction["Interaction States"].items()) {
                    SimulationParameters::InteractionStateParameters interaction_state_parameters{};
                    interaction_state_parameters.name = state.key();
                    interaction_state_parameters.interaction_type = state.value().value("type", "InteractionType");
                    interaction_state_parameters.must_overhead = state.value().value("must_overhead", 0.0);
                    interaction_state_parameters.adhere = state.value().value("adhere", false);

                    if (state.value().find("next_states") != state.value().end()) {
                        for (const auto &next_state:state.value()["next_states"].items()) {
                            interaction_state_parameters.next_states.emplace_back(next_state.key(),
                                                                                  next_state.value()["rate"]);
                        }
                    }
                    interaction_parameters->states.emplace_back(interaction_state_parameters);
                }
                if (auto conditions = interaction.find("Conditions");conditions != interaction.end()) {
                    for (const auto&[key, values]:conditions->items()) {
                        interaction_parameters->cell_conditions.emplace_back(key,
                                                                             values.get<std::vector<std::string>>());
                    }
                }
                parameters.interaction_parameters.emplace_back(std::move(interaction_parameters));
            }
        }
        json_file.close();
        return parameters;
    }

    AnalyserParameters getAnalyserParameters(const std::string &analyser_config) {
        AnalyserParameters parameters{};
        std::ifstream json_file(analyser_config);
        json json_parameters;
        json_file >> json_parameters;
        parameters.active_measurements = json_parameters["Agent-Based-Framework"]["active_measurements"].get<std::unordered_set<std::string>>();
        parameters.cell_state_count = json_parameters["Agent-Based-Framework"]["cell_state_count"].get<std::vector<std::string>>();
        json_file.close();
        return parameters;
    }

    VisualizerParameters getViualizerParameters(const std::string &visualizer_config) {
        VisualizerParameters parameters{};
        std::ifstream json_file(visualizer_config);
        json json_parameters;
        json_file >> json_parameters;

        parameters.pov_active = json_parameters["Agent-Based-Framework"].value("activated", false);
        parameters.run_id = json_parameters["Agent-Based-Framework"].value("runId", 1);
        parameters.output_interval = json_parameters["Agent-Based-Framework"].value("output_interval", 100);
        parameters.include_time = json_parameters["Agent-Based-Framework"].value("includeTimestamp", true);
        parameters.output_video = json_parameters["Agent-Based-Framework"].value("outputVideo", false);
        parameters.px_width = json_parameters["Agent-Based-Framework"].value("pxWidth", "2000");
        parameters.px_height = json_parameters["Agent-Based-Framework"].value("pxHeight", "1900");
        parameters.path_style_transfer = json_parameters["Agent-Based-Framework"].value("path_style_transfer", "");
        parameters.white_noise = json_parameters["Agent-Based-Framework"].value("white_noise", 0.0);
        parameters.image_noise = json_parameters["Agent-Based-Framework"].value("ImageNoise", false);
        parameters.camera_angle = json_parameters["Agent-Based-Framework"].value("camera-angle", 0.0);
        auto position = json_parameters["Agent-Based-Framework"]["camera-position"].get<std::vector<double>>();
        parameters.camera_position = {position[0], position[1], position[2]};
        auto look_at = json_parameters["Agent-Based-Framework"]["camera-look_at"].get<std::vector<double>>();
        parameters.camera_look_at = {look_at[0], look_at[1], look_at[2]};
        for (auto light:json_parameters["Agent-Based-Framework"]["Lightsources"]) {
            parameters.light_sources.emplace_back(Coordinate3D{light[0], light[1], light[2]});
        }
        json_file.close();
        return parameters;
    }

    InputParameters getInputParameters(const std::string &input_config) {
        InputParameters parameters{};
        std::ifstream json_file(input_config);
        json json_parameters;
        json_file >> json_parameters;
        for (auto &distribution:json_parameters["Agent-Based-Framework"]["Distributions"].items()) {
            InputParameters::DistributionParameters para_dist{};
            para_dist.key = distribution.key();
            para_dist.type = distribution.value()["type"];
            para_dist.lambda = distribution.value().value("lambda", 0.0);
            para_dist.source_file = distribution.value().value("source_file", "");
            para_dist.number_of_distributions = distribution.value().value("number-of-distributions", 0);
            parameters.distributions.emplace_back(para_dist);
        }
        for (auto &rate:json_parameters["Agent-Based-Framework"]["Rates"].items()) {
            std::unique_ptr<InputParameters::DefaultRateParameters> para_rate;
            const auto type = rate.value().value("type", "ConstantRate");
            para_rate = std::make_unique<InputParameters::DefaultRateParameters>();
            para_rate->type = type;
            para_rate->key = rate.key();
            para_rate->rate = rate.value().value("rate", 0.0);
            para_rate->condition = rate.value().value("condition", "");
            parameters.rates.emplace_back(std::move(para_rate));
        }
        json_file.close();
        return parameters;
    }

    void executeShellCommand(const std::string &command, const bool suppress_output) {
        if (auto system_return = std::system(command.c_str());!suppress_output) {
            DEBUG_STDOUT("shell.command: " + command + "with return " + std::to_string(system_return));
        }
    }

}