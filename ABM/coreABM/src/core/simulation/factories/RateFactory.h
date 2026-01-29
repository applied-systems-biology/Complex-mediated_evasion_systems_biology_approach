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

#ifndef CORE_SIMULATION_RATEFACTORY_H
#define CORE_SIMULATION_RATEFACTORY_H

#include <boost/filesystem.hpp>

#include "core/utils/io_util.h"

class Rate;
class RateFactory {
    using RateParameter = std::unique_ptr<abm::util::InputParameters::DefaultRateParameters>;
public:
    /// Factory class for rate related calculations defining interaction rates. Rate definitions are given in the associated input-config.
    RateFactory(const std::vector<RateParameter> &rates);

    Rate *getRate(const std::string &rate_key);
    std::unordered_map<std::string, std::unique_ptr<Rate>> &getRates();
protected:
    std::unordered_map<std::string, std::unique_ptr<Rate>> rates_;
};

#endif /* CORE_SIMULATION_RATEFACTORY_H */
