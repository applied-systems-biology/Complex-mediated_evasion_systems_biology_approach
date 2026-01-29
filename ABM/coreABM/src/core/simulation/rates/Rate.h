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

#ifndef CORE_SIMULATION_RATE_H
#define CORE_SIMULATION_RATE_H

#include <string>

class Condition;
class Cell;
class Site;

class Rate {
  public:
    /// Class for rate related calculations defining interaction rates. Rate definitions are given in the associated input-config.
    virtual ~Rate() = default;

    [[nodiscard]] virtual double
    calculateProbability(double timestep, Condition *cond, Cell *cell, Site *site) = 0;
    [[nodiscard]] virtual double getRateValue() const = 0;
    [[nodiscard]] virtual std::string_view getRateType() const = 0;
    virtual void adjustRate(const std::string &, double timestep_size) {};
};

#endif /* CORE_SIMULATION_RATE_H */
