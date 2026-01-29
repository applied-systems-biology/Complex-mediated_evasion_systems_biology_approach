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

#ifndef RATEFACTORYCME_H
#define RATEFACTORYCME_H
#include <core/simulation/factories/RateFactory.h>


class RateFactoryCME: public RateFactory {
    using RateParameter = std::unique_ptr<abm::util::InputParameters::DefaultRateParameters>;
public:
    RateFactoryCME(const std::vector<RateParameter> &rates, const double timestep);
};



#endif //RATEFACTORYCME_H
