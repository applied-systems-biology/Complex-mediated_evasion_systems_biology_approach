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

#include "core/visualisation/PovFile.h"
#include "core/utils/io_util.h"
#include "core/simulation/Site.h"
#include "core/visualisation/PovRayObject.h"
#include "core/utils/macros.h"

PovFile::PovFile(const std::string &filebody) {
    filebody_pov_ = filebody;
}

void PovFile::setGlobalPart(std::string global) {
    global_part_ = global;
}

void PovFile::setCameraPart(std::string cam) {
    camera_part_ = cam;
}

void PovFile::setBackgroundPart(std::string bg) {
    background_part_ = bg;
}

void PovFile::setLightPart(std::string light) {
    light_part_ = light;
}

void PovFile::addLightPart(const std::string &light) {
    light_part_ += light;
}

void PovFile::setBorderPart(std::string border) {
    border_part_ = border;
}

void PovFile::addPovObject(const std::string &povObj) {
    pov_objects_.emplace_back(povObj);
}

void PovFile::setDimensions(const std::string px_width, const std::string px_height) {
    px_width_ = px_width, px_height_ = px_height;
}

void PovFile::setIncludeTime(bool include_time) {
    include_time_ = include_time;
}

void PovFile::doPovProcess(double current_time, std::string zpos) {
    std::ofstream output(filebody_pov_ + zpos + ".pov");
    output << global_part_ << '\n';
    output << camera_part_ << '\n';
    output << light_part_ << '\n';
    output << border_part_ << '\n';
    output << background_part_ << '\n';
    for (const auto &obj: pov_objects_) {
        output << obj << '\n';
    }
    output.close();
    int left_padding = static_cast<int>(stoi(px_width_) / 5);
    int below_padding = 50;
    int font_size = static_cast<int>(stoi(px_height_) / 30);
    std::string font_color = "black";
    int X = std::stoi(px_width_) - left_padding;
    int Y = std::stoi(px_height_) - below_padding;
    char timeOfFrame[20];
    sprintf(timeOfFrame, "%.3f", current_time);
    std::ostringstream command;
    command << "povray -GA -d +H" << px_height_ << " +W" << px_width_ << " +I" << filebody_pov_ << zpos << ".pov +O"
            << filebody_pov_ << zpos << ".png >/dev/null 2>&1 ";
    if (include_time_) {
        command << "&& convert -pointsize " << font_size << " -fill " << font_color << " -draw \"text " << X << " " << Y
                << " '" << timeOfFrame
                << " s'\" " << filebody_pov_ << zpos << ".png " << filebody_pov_ << ".png ";
    }
    command << "&& rm -f " << filebody_pov_ << zpos << ".pov";

    abm::util::executeShellCommand(command.str(), true);
}

void PovFile::transcribeSite(Site &site) {
    // set background color
    addPovObject("background { color rgb <0.9, 0.9, 0.9> }");

    if (abm::util::isSubstring("CuboidSite", site.getType())) {
        auto lowerBoundReal = site.getLowerLimits();
        auto upperBoundReal = site.getUpperLimits();
        double extension = 0.1;
        double thickness = 0.05;
        double threshold_2d = 2.0;
        Coordinate3D bound_extension = {1.0, 1.0, 1.0};
        bound_extension *= extension;
        auto lowerBound = lowerBoundReal - bound_extension;
        auto upperBound = upperBoundReal + bound_extension;
        Coordinate3D point1, point2, point3, point4, point5, point6, point7, point8;
        point1 = lowerBound;
        point2 = Coordinate3D{upperBound.x, lowerBound.y, lowerBound.z};
        point3 = Coordinate3D{upperBound.x, upperBound.y, lowerBound.z};
        point4 = Coordinate3D{lowerBound.x, upperBound.y, lowerBound.z};
        point5 = Coordinate3D{lowerBound.x, upperBound.y, upperBound.z};
        point6 = Coordinate3D{lowerBound.x, lowerBound.y, upperBound.z};
        point7 = Coordinate3D{upperBound.x, lowerBound.y, upperBound.z};
        point8 = upperBound;

        ColorRGB colorBounds = ColorRGB(0.3, 0.3, 0.3, 0.3);//(0.5, 0.5, 0.5, 0.5);
        if (abs(lowerBoundReal.z - upperBoundReal.z) < threshold_2d) {
            //two-dimensional case
            addPovObject(PovRayObject::getCylinder(point1, point2, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point2, point3, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point3, point4, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point1, point4, thickness, colorBounds));
        } else {
            //three-dimensional case
            addPovObject(PovRayObject::getCylinder(point1, point2, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point2, point3, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point3, point4, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point1, point4, thickness, colorBounds));

            addPovObject(PovRayObject::getCylinder(point5, point6, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point6, point7, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point7, point8, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point8, point5, thickness, colorBounds));

            addPovObject(PovRayObject::getCylinder(point4, point5, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point1, point6, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point2, point7, thickness, colorBounds));
            addPovObject(PovRayObject::getCylinder(point3, point8, thickness, colorBounds));
        }
    }
}

void PovFile::transcribeAgents(Site &site) {

    for (auto agent: *site.getAgentManager()) {
        for (auto sphere: agent->getSurface()->getAllSpheresOfThis()) {
            auto position = sphere->getPosition();
            auto radius = sphere->getRadius();
            auto color = agent->getSurface()->getColorRGB();
            addPovObject(PovRayObject::getSphere(position, radius, *color));
        }
    }
}