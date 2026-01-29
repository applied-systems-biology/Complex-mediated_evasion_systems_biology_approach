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

#ifndef CORE_SIMULATION_CuboidSiteCME_H
#define CORE_SIMULATION_CuboidSiteCME_H

#include <string>
#include <iostream>

#include "InSituMeasurementsCME.h"
#include "core/simulation/Site.h"
#include "io_utilsCME.h"

class CuboidSiteCME : public Site {
 public:

  CuboidSiteCME(Randomizer *random_generator,
                    std::shared_ptr<InSituMeasurements> measurements,
             std::string config_path, std::unordered_map<std::string,
             std::string> cmd_input_args, std::string output_path);
  ~CuboidSiteCME() override = default;
  void handleBoundaryCross(Agent * agent, Coordinate3D * movement, double current_time) final;
  [[nodiscard]] bool containsPosition(Coordinate3D) final;
  [[nodiscard]] std::string getType() const override{return "CuboidSiteCME";}
  Coordinate3D getRandomPosition(double cell_radius);
  Coordinate3D getRandomPosition() final {};

  Coordinate3D getLowerLimits() final;
  Coordinate3D getUpperLimits() final;
  Coordinate3D getRandomBoundaryPoint() final;

  Coordinate3D generateRandomDirectionVector(Coordinate3D position, double length) final;
  std::vector<Coordinate3D> getSystemBoundaries() {return {lower_bound_, upper_bound_};};
  Coordinate3D generatePersistentDirectionVector(Coordinate3D position,
                                                 double length,
                                                 Coordinate3D prev_vector,
                                                 double previous_alpha) final;

  Coordinate3D generateBackShiftOnContacting(SphereRepresentation *active_sphere,
                                             SphereRepresentation *passiveSphere,
                                             double must_overhead) override;

  void initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                                      double current_time,
                                      double time_delta);

  Coordinate3D random_point_on_sphere(double radius, Coordinate3D pos);
  Coordinate3D biased_random_point_on_sphere(double radius, Coordinate3D pos_of_cell, Coordinate3D pos_of_uptake);

  void secretion_of_agents_at_cell_surface_randomly(double timestep, double current_time, double rate, const std::string& agent_type, Coordinate3D pos, double radius);
  void secretion_of_agents_at_cell_surface(double timestep, double current_time, double rate, const std::string& agent_type, Coordinate3D pos_of_cell, double radius, std::vector<Coordinate3D> pos_of_uptake);

  void doAgentDynamics(Randomizer *random_generator, SimulationTime &time) override;

  Coordinate3D getCenterPosition();

  void handleCmdInputArgs(std::unordered_map<std::string, std::string> cmd_input_args) override;


  void do_site_dynamics(double timestep, double current_time) override;

  Coordinate3D two_cells_layout(double radius, int cell, double pos_dist);

  bool check_steady_state(double dt);

  bool checkForStopping() const override final {return stop_sim;}

protected:
  Coordinate3D upper_bound_;
  Coordinate3D lower_bound_;
  std::tuple<int, int, int> molecules_grid_size_;
  Coordinate3D getZEqualsZeroPosition(double radius);
  bool check_collision(Coordinate3D pos_to_check);

 std::map<std::string, int> prev_populations;

  std::map<std::string, double> flow_;

  abm::utilCME::StoppingCriteria stopping_criteria_;

  bool stop_sim{};
};

#endif /* CORE_SIMULATION_CuboidSiteCME_H */