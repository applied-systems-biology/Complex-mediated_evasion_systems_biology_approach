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

#ifndef CORE_SIMULATION_RANDOMWALK_H
#define CORE_SIMULATION_RANDOMWALK_H

#include "Movement.h"

#include "core/basic/Randomizer.h"
#include "core/basic/InversionSampler.h"

class RandomWalk : public Movement {
public:
  // Class for the normal random walk
    RandomWalk(Site *site, unsigned int spatial_dimensions);
    RandomWalk(Site *site, double speed, unsigned int spatial_dimensions);
    RandomWalk(Site *site, double speed, double speedStddev, unsigned int spatial_dimensions);

    Coordinate3D *move(double timestep, double dc);
    Coordinate3D *calculateDiffusiveMove(double, double dc);
    Coordinate3D *calculateRandomMove(double);
    double getSpeed() { return vector_length_per_timeunit_; };
    void setSpeed(double speed) { vector_length_per_timeunit_ = speed; };
    std::string getMovementName() final { return "RandomWalk"; }

private:
    double vector_length_per_timeunit_{};
    double speed_stddev_{};
};

#endif /* CORE_SIMULATION_RANDOMWALK_H */
