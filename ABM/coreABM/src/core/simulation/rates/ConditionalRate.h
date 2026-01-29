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

#ifndef CORE_SIMULATION_CONDITIONALRATE_H
#define CORE_SIMULATION_CONDITIONALRATE_H

#include <memory>

#include "Rate.h"
#include "core/simulation/Condition.h"

class ConditionalRate : public Rate {
public:
  // Class for a constant rate that is activated for a specific conditioned specified in the input configuration
    ConditionalRate(double constant_value, std::unique_ptr<Condition> condition) : condition_(std::move(condition)),
                                                                                   constant_rate_(constant_value) {}

    [[nodiscard]] double calculateProbability(double timestep, Condition *cond, Cell *cell, Site *site) final;
    [[nodiscard]] std::string_view getRateType() const final { return "ConditionalRate"; }
    [[nodiscard]] double getRateValue() const final { return constant_rate_; }

private:
    std::unique_ptr<Condition> condition_;
    double constant_rate_;
};

#endif /* CORE_SIMULATION_CONDITIONALRATE_H */
