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

#include <sstream>

#include "core/basic/ColorRGB.h"

ColorRGB::ColorRGB() {
    red = 0;
    green = 0;
    blue = 0;

    transmit = 0;
    this->colorName = "";
}

ColorRGB::ColorRGB(double redPart, double greenPart, double bluePart) {
    if (redPart <= 1 && redPart >= 0) {
        red = redPart;
    } else {
        red = 0;
    }

    if (greenPart <= 1 && greenPart >= 0) {
        green = greenPart;
    } else {
        green = 0;
    }

    if (bluePart <= 1 && bluePart >= 0) {
        blue = bluePart;
    } else {
        green = 0;
    }
    this->colorName = "";
    transmit = 0;
}

ColorRGB::ColorRGB(double redPart, double greenPart, double bluePart, double transmit) {
    if (redPart <= 1 && redPart >= 0) {
        red = redPart;
    } else {
        red = 0;
    }

    if (greenPart <= 1 && greenPart >= 0) {
        green = greenPart;
    } else {
        green = 0;
    }

    if (bluePart <= 1 && bluePart >= 0) {
        blue = bluePart;
    } else {
        green = 0;
    }
    this->colorName = "";
    this->transmit = transmit;
}

ColorRGB::ColorRGB(Randomizer *randomizer, std::string colorName) {
    red = 0;
    green = 0;
    blue = 0;
    transmit = 0;
    if (colorName == "random") {
        setRandomColor(randomizer);
    } else {
        setColor(randomizer, colorName);
    }
    this->colorName = colorName;
}

std::string ColorRGB::printPovColorRGB() {

    std::ostringstream ss;
    ss << "rgb <" << red << "," << green << "," << blue << ">";

    return ss.str();
}

std::string ColorRGB::printPovColorRGBF() {

    std::ostringstream ss;
    ss << "rgbf <" << red << "," << green << "," << blue << ",0>";

    return ss.str();
}

std::string ColorRGB::printPovColorRGBT() {
    std::ostringstream ss;
    ss << "rgbt <" << red << "," << green << "," << blue << "," << transmit << ">";

    return ss.str();
}

void ColorRGB::setColor(Randomizer *randomizer, std::string colorName) {

    if (!colorName.compare("red")) red = 1;
    if (!colorName.compare("green")) green = 1;
    if (!colorName.compare("blue")) blue = 1;
    if (!colorName.compare("white")) {
        red = 1;
        green = 1;
        blue = 1;
    }
    if (!colorName.compare("black")) {
        red = 0;
        green = 0;
        blue = 0;
    }
    if (!colorName.compare("orange")) {
        red = 1;
        green = 0.4;
    }
    if (!colorName.compare("pink")) {
        red = 1;
        green = 0.08;
        blue = 0.58;
    }
    if (!colorName.compare("violet")) {
        red = 0.58;
        green = 0.0;
        blue = 0.83;
    }
    if (!colorName.compare("cyan")) {
        red = 0;
        green = 1;
        blue = 1;
    }
    if (!colorName.compare("aquamarin")) {
        red = 0.49;
        green = 1;
        blue = 0.83;
    }
    if (!colorName.compare("gold")) {
        red = 1;
        green = 0.84;
        blue = 0;
    }
    if (!colorName.compare("grey")) {
        red = 0.5;
        green = 0.5;
        blue = 0.5;
    }
    if (!colorName.compare("greenTransp")) {
        //green=rand.generateDouble(); red=rand.generateDouble(); blue = rand.generateDouble(); transmit=0.35 + rand.generateDouble(0.25);
        green = 1;
        transmit = 0.35 + randomizer->generateDouble(0.25);
    }
    if (!colorName.compare("darkGreen")) {
        //green=rand.generateDouble(); red=rand.generateDouble(); blue = rand.generateDouble(); transmit=0.35 + rand.generateDouble(0.25);
        green = 0.4;
        transmit = 0.5 + randomizer->generateDouble(0.25);
    }
    if (!colorName.compare("redTransp")) {
        red = 1;
        transmit = 0.5 + randomizer->generateDouble(0.25);
    }
    if (!colorName.compare("blueTransp")) {
        blue = 1;
        transmit = 0.35 + randomizer->generateDouble(0.25);
    }
    if (!colorName.compare("blueTransp2")) {
        red = 0.12;
        green = 0.56;
        blue = 1;
        transmit = 0.45 + randomizer->generateDouble(0.25);
    }
    if (!colorName.compare("orangeTransp")) {
        red = 1;
        green = 0.4;
        transmit = 0.35 + randomizer->generateDouble(0.25);
    }

    return;

}

void ColorRGB::setRandomColor(Randomizer *randomizer){
    this->red = randomizer->generateDouble(1.0);
    this->green = randomizer->generateDouble(1.0);
    this->blue = randomizer->generateDouble(1.0);
    this->transmit = 0.5 + randomizer->generateDouble(0.5);
}

void ColorRGB::setColor(double red, double green, double blue) {
    this->red = red;
    this->green = green;
    this->blue = blue;
}

void ColorRGB::setColor(double red, double green, double blue, double transparency) {
    this->red = red;
    this->green = green;
    this->blue = blue;
    this->transmit = transparency;
}

double ColorRGB::getRed() const {
    return red;
}

double ColorRGB::getGreen() const {
    return green;
}

double ColorRGB::getBlue() const {
    return blue;
}

double ColorRGB::getTransmit() const {
    return transmit;
}

std::string ColorRGB::getColorName() const {
    return colorName;
}