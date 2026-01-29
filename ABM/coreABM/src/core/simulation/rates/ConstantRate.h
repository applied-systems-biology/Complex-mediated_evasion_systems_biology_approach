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

#ifndef CORE_SIMULATION_CONSTANTRATE_H
#define CORE_SIMULATION_CONSTANTRATE_H

#include "Rate.h"

class ConstantRate : public Rate {
public:
  // Class for a constant rate which does not change during runtime.
    explicit ConstantRate(double constant_value) : constant_rate_(constant_value) {}

    [[nodiscard]] double calculateProbability(double timestep, Condition *cond, Cell *cell, Site *site) final;
    [[nodiscard]] double getRateValue() const final { return constant_rate_; }
    [[nodiscard]] std::string_view getRateType() const final { return "ConstantRate"; }

private:
    double constant_rate_;
};
#endif /* CORE_SIMULATION_CONSTANTRATE_H */
