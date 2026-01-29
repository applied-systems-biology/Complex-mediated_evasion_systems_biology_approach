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

#include "core/simulation/Condition.h"
#include "core/simulation/Cell.h"

Cell *Condition::getCell() {
    return cell_;
}

std::string Condition::getStringCondition() {
    return condition_;
}

bool Condition::isFulfilled(Condition *condition) {
    bool fulfilled = false;
    if (cell_ != nullptr) {
        if (condition->getCell() == nullptr) {
            fulfilled = (condition->getStringCondition() == cell_->getTypeName());
        } else {
            fulfilled = (cell_ == condition->getCell());
        }

    } else {
        if (condition->getCell() == nullptr) {
            fulfilled = (condition->getStringCondition() == condition_);
        } else {
            fulfilled = (condition->getCell()->getTypeName() == condition_);
        }
    }
    return fulfilled;
}
