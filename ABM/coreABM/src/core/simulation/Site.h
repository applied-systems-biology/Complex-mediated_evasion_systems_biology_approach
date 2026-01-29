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

#ifndef CORE_SIMULATION_SITE_H
#define CORE_SIMULATION_SITE_H

#include <memory>
#include <utility>
#include <vector>
#include <algorithm>

#include "core/utils/time_util.h"
#include "core/basic/Coordinate3D.h"
#include "core/simulation/Cell.h"
#include "core/analyser/InSituMeasurements.h"
#include "core/simulation/morphology/SphereRepresentation.h"
#include "core/simulation/neighbourhood/BoundaryCondition.h"
#include "core/simulation/neighbourhood/NeighbourhoodLocator.h"
#include "core/simulation/AgentManager.h"
#include "core/simulation/factories/RateFactory.h"
#include "core/simulation/factories/CellFactory.h"
#include "core/simulation/factories/InteractionFactory.h"
#include "core/simulation/factories/CellStateFactory.h"
#include "core/simulation/factories/InteractionStateFactory.h"


class Agent; //forward declaration
class Analyser;
class Randomizer;
class CellFactory;
class InteractionStateFactory;

class Site {

public:
    /// Class for handling environment interactions and wrapping functionality of all aspects that are happening during the simulation inside of the environment (i.e. site)
    Site(Randomizer *random_generator,
         std::shared_ptr<InSituMeasurements> measurements,
         std::string config_path,
         std::unordered_map<std::string, std::string> cmd_input_args,
         std::string output_path);

    virtual ~Site() = default;

    /*!
     * Conducts all actions that happen in one timestep in the site
     * @param random_generator Randomizer object that contains randomizer of current run
     * @param time SimulationTime object, i.e. contains current time and timestep
     */
    virtual void doAgentDynamics(Randomizer *random_generator, SimulationTime &time);

    /*!
     * Sets stopping conditions of simulator-config.json
     * @param stopping_criteria Vector of Strings of stopping criteria (i.e. "stopping_criteria": ["FirstPassageTime"])
     */
    //void setStoppingCondition(const std::vector<std::string> &stopping_criteria);

    /*!
     * Terminates Simulation for certain interactions (i.e. all phagocytes were touched)
     * @param interaction Interaction object that contains an interaction
     */
    //void stopRunForCertainState(Cell &cell, std::string state, double current_time);

    virtual void handleCmdInputArgs(std::unordered_map<std::string, std::string> cmd_input_args);
    void addOutputPath(const std::string& path) {visualizer_output_paths_.emplace_back(path);};
    void addOutputCommand(const std::string& path) {output_commands_.emplace_back(path);}
    virtual void receiveFrontendParameter(abm::util::SimulationParameters& sim_para, abm::util::InputParameters& inp_para,
                                  const std::string &config_path, const std::string &output_path, std::string sid);
    std::vector<std::string> getOutputPaths() {return visualizer_output_paths_;};
    std::vector<std::string> getOutputCommands() {return output_commands_;};
    virtual bool checkForStopping() const;
    Randomizer *getRandomGenerator() { return random_generator_; }
    NeighbourhoodLocator *getNeighbourhoodLocator() { return neighbourhood_locator_.get(); }
    AgentManager *getAgentManager() const { return agent_manager_.get(); }
    InSituMeasurements *getMeasurments() const { return measurements_.get(); }
    RateFactory *getRateFactory() const { return rate_factory_.get(); }
    CellFactory *getCellFactory() const { return cell_factory_.get(); }
    InteractionFactory *getInteractionFactory() const { return interaction_factory_.get(); }
    CellStateFactory *getCellStateFactory() const { return cell_state_factory_.get(); }
    InteractionStateFactory *getInteractionStateFactory() const { return interaction_state_factory_.get(); }
    Coordinate3D getBoundaryInputVector() { return boundary_input_vector_; }
    bool getLargeTimestepActive() { return large_timestep_active; }
    double getMaxTime() {return parameters_.max_time;}
    double getTimeStepping() {return parameters_.time_stepping;}
    abm::util::VisualizerParameters getOverwrittenVisParameters() {return parameters_.visualizer_to_overwrite;};
    [[nodiscard]] int getState() const { return state_; }
    [[nodiscard]] unsigned int getNumberOfSpatialDimensions() const { return dimensions; }
    [[nodiscard]] double getLatestAlpha2dTurningAngle() const { return alpha2dTurningAngle; }
    [[nodiscard]] double getInputRate(std::string agent_name) { return input_rates_[agent_name]; }
    [[nodiscard]] std::string getIdentifier() const { return identifier_; }

    friend void InSituMeasurements::observeMeasurements(const SimulationTime &time);

    virtual void handleBoundaryCross(Agent *, Coordinate3D *, double current_time) = 0;
    virtual bool containsPosition(Coordinate3D position) = 0;
    [[nodiscard]] virtual std::string getType() const = 0;
    virtual Coordinate3D getRandomPosition() = 0;
    virtual Coordinate3D getRandomBoundaryPoint() = 0;
    virtual Coordinate3D getLowerLimits() = 0;
    virtual Coordinate3D getUpperLimits() = 0;
    virtual Coordinate3D generateRandomDirectionVector(Coordinate3D position, double length) = 0;
    virtual Coordinate3D generatePersistentDirectionVector(Coordinate3D position,
                                                           double length,
                                                           Coordinate3D prevVector,
                                                           double previousAlpha) = 0;
    virtual Coordinate3D generateBackShiftOnContacting(SphereRepresentation *activeSphere,
                                                       SphereRepresentation *passiveSphere,
                                                       double mustOverhead) = 0;
    virtual double getFeatureValueByName(std::string name) { return {}; }
    [[nodiscard]] virtual double getRadius() const { return 0.0; }
    virtual Coordinate3D generateDirectedVector(Coordinate3D position,
                                                SphericCoordinate3D posOfGoal,
                                                double length) { return {}; }
    virtual Coordinate3D generateDirectedVector(Coordinate3D position,
                                                double alpha, double length) { return {}; }
    virtual std::vector<Coordinate3D> getSystemBoundaries() { return {}; };


    virtual void do_site_dynamics(double timestep, double current_time){};

    std::string get_boundary_condition(){return boundary_type_;}

    BoundaryCondition* get_boundary(){return boundary_condition_.get();}

  protected:
    void setBoundaryCondition();
    virtual void initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                          const std::string &input_dir,
                          double current_time,
                          double time_delta){};
    int state_{};
    unsigned int dimensions{};
    std::map<std::string, double> input_rates_{};
    double alpha2dTurningAngle{};
    std::string identifier_{};
    std::string boundary_type_{};
    bool large_timestep_active = false;
    std::vector<int> detected_fungi_id{};
    Coordinate3D boundary_input_vector_{};
    std::vector<std::pair<std::string, long>> stopping_cell_states;
    Randomizer *random_generator_;
    std::shared_ptr<InSituMeasurements> measurements_;
    std::unique_ptr<BoundaryCondition> boundary_condition_;
    std::unique_ptr<NeighbourhoodLocator> neighbourhood_locator_;
    std::unique_ptr<AgentManager> agent_manager_;

    std::string config_path_{};
    std::vector<std::string> visualizer_output_paths_{};
    std::vector<std::string> output_commands_{};
    abm::util::SimulationParameters parameters_{};
    abm::util::InputParameters input_parameters_{};

    std::unique_ptr<RateFactory> rate_factory_;
    std::unique_ptr<CellFactory> cell_factory_;
    std::unique_ptr<InteractionFactory> interaction_factory_;
    std::unique_ptr<CellStateFactory> cell_state_factory_;
    std::unique_ptr<InteractionStateFactory> interaction_state_factory_;
};

#endif    /* CORE_SIMULATION_SITE_H */

