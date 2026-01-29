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

#ifndef CORE_VISUALISATION_POVFILE_H
#define CORE_VISUALISATION_POVFILE_H

#include <string>
#include <vector>

class Site;

class PovFile {
public:
    explicit PovFile(const std::string &filebody);

    void setGlobalPart(std::string);
    void setCameraPart(std::string);
    void setLightPart(std::string);
    void addLightPart(const std::string &);
    void setDimensions(std::string px_width, std::string px_height);
    void setBorderPart(std::string);
    void setBackgroundPart(std::string);
    void addPovObject(const std::string &);
    void setIncludeTime(bool);
    void doPovProcess(double current_time, std::string zpos = "");
    virtual void transcribeSite(Site &site);
    virtual void transcribeAgents(Site &site);

private:
    std::string global_part_;
    std::string camera_part_;
    std::string light_part_;
    std::string px_width_;
    std::string px_height_;
    std::string border_part_;
    std::string background_part_;
    std::string filebody_pov_;
    std::vector<std::string> pov_objects_;
    bool include_time_;
};

#endif /* CORE_VISUALISATION_POVFILE_H */
