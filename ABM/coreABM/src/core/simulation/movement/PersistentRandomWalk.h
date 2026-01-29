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

#ifndef CORE_SIMULATION_PERSISTENTRANDOMWALK_H
#define CORE_SIMULATION_PERSISTENTRANDOMWALK_H

#include <cmath>

#include "Movement.h"
#include "core/basic/Randomizer.h"
#include "core/basic/Sampler.h"


class Agent;

class PersistentRandomWalk : public Movement {
public:
    PersistentRandomWalk(Agent *, double, double, unsigned int spatial_dimensions);

    Coordinate3D *move(double, double diffusionConstant);
    double getStartingTime() final;
    void setPreviousMove(Coordinate3D *) final;
    void annulatePersistence() final{ persistence_time_left_ = 0; }

    Coordinate3D *movePersistent(double);
    Coordinate3D *moveRandomly(double);
    void decrementLeftTime(double);
    void setNewPersistence();
    [[nodiscard]] bool persistentMove() const;
    void setSpeed(double speed) { speed_ = speed; };

protected:
    Agent *agent_{};
    double speed_{};
    std::unique_ptr<Coordinate3D> current_velocity_{};
    double persistent_angle_alpha_2_d_{};

private:
    double persistence_time_{};
    double persistence_time_left_{};
    double persistence_time_start_{};
    std::unique_ptr<Coordinate3D> persistence_direction_{};
    std::shared_ptr<Sampler> sampler_{};

};

#endif /* CORE_SIMULATION_PERSISTENTRANDOMWALK_H */
