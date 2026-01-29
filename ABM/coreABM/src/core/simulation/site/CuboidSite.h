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

#ifndef CORE_SIMULATION_CUBOIDSITE_H
#define CORE_SIMULATION_CUBOIDSITE_H

#include <string>
#include <iostream>

#include "core/simulation/Site.h"


class CuboidSite : public Site {
 public:

  CuboidSite(Randomizer *random_generator,
             std::shared_ptr<InSituMeasurements> measurements,
             std::string config_path, std::unordered_map<std::string,
             std::string> cmd_input_args, std::string output_path);

  void handleBoundaryCross(Agent * agent, Coordinate3D * movement, double current_time) final;
  [[nodiscard]] bool containsPosition(Coordinate3D) final;
  [[nodiscard]] std::string getType() const override{return "CuboidSite";}
  Coordinate3D getRandomPosition() final;
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
 protected:
  Coordinate3D upper_bound_;
  Coordinate3D lower_bound_;
  std::tuple<int, int, int> molecules_grid_size_;
};

#endif /* CORE_SIMULATION_CUBOIDSITE_H */