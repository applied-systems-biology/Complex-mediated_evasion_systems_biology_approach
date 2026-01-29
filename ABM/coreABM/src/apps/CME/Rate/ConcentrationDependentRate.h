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

//
// Created by ybachelot on 15.01.25.
//

#ifndef CONCENTRATIONDEPENDENTRATE_H
#define CONCENTRATIONDEPENDENTRATE_H


#include <core/utils/io_util.h>

#include "core/simulation/rates/Rate.h"

class ConcentrationDependentRate: public Rate {
public:
    explicit ConcentrationDependentRate(const double initial_rate, const double alpha) : initial_rate_(initial_rate), alpha_(alpha), current_rate_(initial_rate){};

    [[nodiscard]] double calculateProbability(double timestep, Condition *cond, Cell *cell, Site *site) final;
    [[nodiscard]] double getRateValue() const final { return current_rate_; }
    [[nodiscard]] std::string_view getRateType() const final { return "ConcentrationDependentRate"; }

private:
    double initial_rate_{};
    double alpha_{};
    double current_rate_{};
};



#endif //CONCENTRATIONDEPENDENTRATE_H
