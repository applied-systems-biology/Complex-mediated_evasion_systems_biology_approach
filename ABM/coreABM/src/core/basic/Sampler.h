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

#ifndef CORE_BASIC_SAMPLER_H
#define CORE_BASIC_SAMPLER_H

#include <iostream>
#include <vector>

#include "core/basic/Randomizer.h"

class Sampler {
public:
  // Class for randomly sampling values
    Sampler();

    Sampler(const char *);
    virtual ~Sampler() = default;
    virtual void setSampleFunction(int);
    virtual void setSampleRange(double, double, unsigned int dim = 1);
    virtual void sample(Randomizer *randomizer, unsigned int distrColumn = 0);
    virtual double getSampledValue(unsigned int dim = 1);

protected:
    int samplingBasis;
    const char *filename;
    int samplingFunction;
    std::vector<double> lowerBound;
    std::vector<double> upperBound;

};

#endif /* CORE_BASIC_SAMPLER_H */
