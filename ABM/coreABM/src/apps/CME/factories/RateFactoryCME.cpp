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

#include "RateFactoryCME.h"

#include <apps/CME/io_utilsCME.h>

#include "apps/CME/Rate/ConcentrationDependentRate.h"
#include "core/simulation/Condition.h"
#include "core/simulation/rates/ConditionalRate.h"
#include "core/simulation/rates/ConstantRate.h"

RateFactoryCME::RateFactoryCME(const std::vector<RateParameter> &rates, const double timestep): RateFactory(rates) {
    for (const auto &rate: rates) {
        if ("ConstantRate" == rate->type) {
            if (rate->key == "uptake") {
                rates_.emplace(rate->key, std::make_unique<ConstantRate>(rate->rate/timestep));
            }
            else {
                rates_.emplace(rate->key, std::make_unique<ConstantRate>(rate->rate));
            }
        } else if ("ConditionalRate" == rate->type) {
            auto condition = std::make_unique<Condition>(rate->condition);
            rates_.emplace(
                    rate->key, std::make_unique<ConditionalRate>(rate->rate, std::move(condition)));
        }
        else if ("ConcentrationDependentRate" == rate->type) {
            auto *conc_rate = dynamic_cast<abm::utilCME::ConcentrationDependentRateParameters *>(rate.get());
            rates_.emplace(rate->key, std::make_unique<ConcentrationDependentRate>(conc_rate->rate, conc_rate->alpha));
        }
    }
}
