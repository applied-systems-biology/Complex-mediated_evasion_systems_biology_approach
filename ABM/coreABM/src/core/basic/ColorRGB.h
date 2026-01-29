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

#ifndef CORE_BASIC_COLORRGB_H
#define CORE_BASIC_COLORRGB_H

#include <string>

#include "core/basic/Randomizer.h"

class ColorRGB {
public:
  // Class for wrapping color-based information
    ColorRGB();

    ColorRGB(double redPart, double greenPart, double bluePart);
    ColorRGB(double redPart, double greenPart, double bluePart, double transmit);
    ColorRGB(Randomizer *randomizer, std::string colorName);

    std::string printPovColorRGB();

    std::string printPovColorRGBF();

    std::string printPovColorRGBT();

    [[nodiscard]] double getRed() const;

    [[nodiscard]] double getGreen() const;

    [[nodiscard]] double getBlue() const;

    [[nodiscard]] double getTransmit() const;

    void setColor(Randomizer *randomizer, std::string colorName);

    void setRandomColor(Randomizer *randomizer);

    void setColor(double red, double green, double blue);

    void setColor(double red, double green, double blue, double transparancy);

    [[nodiscard]] std::string getColorName() const;

private:

    std::string colorName;

    double red;

    double green;

    double blue;

    double transmit;
};

#endif /* CORE_BASIC_COLORRGB_H */
