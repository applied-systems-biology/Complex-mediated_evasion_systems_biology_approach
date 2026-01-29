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

#include <cassert>

#include <boost/filesystem.hpp>

#include "Visualizer.h"
#include "core/utils/macros.h"
#include "core/utils/time_util.h"
#include "core/simulation/Site.h"
#include "core/visualisation/PovFile.h"
#include "core/visualisation/PovRayObject.h"

Visualizer::Visualizer(const std::string &config_path, const std::string &project_dir, int total_runs) {
    auto const visualizer_config = static_cast<boost::filesystem::path>(config_path).append(
            "visualisation-config.json").string();
    assert(!project_dir.empty());
    parameters_ = abm::util::getViualizerParameters(visualizer_config);

    // Initialize visualization folders for all runs
    if (parameters_.pov_active && total_runs > 0) {

        // if camera view is from above
        if (parameters_.camera_position.x == 0 and parameters_.camera_position.y == 0) {
            double ratio_px_microns = stod(parameters_.px_width) / (2 * tan((M_PI * parameters_.camera_angle * 0.5) / 180) *
                                                                    parameters_.camera_position.z);
            DEBUG_STDOUT("Resolution: " << ratio_px_microns << " pixel per micron ---> " << 1 / ratio_px_microns
                                        << " microns per pixel");
        }

        for (int i = 0; i < total_runs; i++) {
            std::string output_path = "povs_run-" + std::to_string(i + 1);
            visualization_path_.emplace_back(
                    static_cast<boost::filesystem::path>(project_dir).append(output_path).string());
            boost::filesystem::create_directories(visualization_path_[i]);
            DEBUG_STDOUT("Visualisation is activated for run id " << i + 1);
        }


    } else if (parameters_.pov_active && parameters_.run_id == 0) {
        ERROR_STDERR("Visualisation is activated but run id is 0");
    }
}

void Visualizer::visualizeCurrentConfiguration(Site &site, const SimulationTime &time, int run,
                                               bool simulation_end) const {
    bool writeout = (0 == time % static_cast<int>(parameters_.output_interval));

    // Generate visualization via PovFiles
    if (parameters_.pov_active && (writeout || simulation_end)) {
        const auto
                povDirTimestep = static_cast<boost::filesystem::path>(visualization_path_[run - 1]).append(
                abm::util::generateUniformString(time.getCurrentTime()));

        std::string filepath = abm::util::getLastNFolders(povDirTimestep.string(), 4) + ".png";
        site.addOutputPath(filepath); //send filepath to frontend api
        auto newPovFile = std::make_shared<PovFile>(povDirTimestep.string());
        const auto boundary = "";

        std::string globals;
        newPovFile->setGlobalPart(globals);
        newPovFile->setDimensions(parameters_.px_width, parameters_.px_height);
        newPovFile->setCameraPart(PovRayObject::getCamera(parameters_.camera_position, parameters_.camera_look_at,
                                                          parameters_.camera_angle));
        newPovFile->setBackgroundPart(PovRayObject::getBackground(ColorRGB(1, 1, 1)));
        newPovFile->setLightPart("");
        for (const auto &lightsource : parameters_.light_sources) {
            newPovFile->addLightPart(PovRayObject::getLightsource(lightsource));
        }
        newPovFile->setIncludeTime(parameters_.include_time);
        newPovFile->setBorderPart(boundary);
        newPovFile->transcribeSite(site);
        newPovFile->transcribeAgents(site);
        newPovFile->doPovProcess(time.getCurrentTime());

        // Conclude simulation and generate video
        if (simulation_end) {
            concludeRun(site, run);
        }
    }
}

void Visualizer::concludeRun(Site &site, int run) const {
    using boost::filesystem::path;
    std::string w = parameters_.px_width, h = parameters_.px_height;

    // Create video of simulated data
    std::ostringstream command;
    command << "ffmpeg -pattern_type glob -i " << boost::filesystem::absolute(static_cast<path>(visualization_path_[run-1])).append("*.png")
            << " -an -vb 6k -s " << w << "x" << h << " -vbt 1M -pass 1 -vcodec png -f mp4 -y "//simulation.mp4";
            << boost::filesystem::absolute(static_cast<path>(visualization_path_[run-1])).append("simulation.mp4")
            << " && yes\"N\" | ffmpeg -i " << boost::filesystem::absolute(static_cast<path>(visualization_path_[run-1])).append("simulation.mp4")
            << " -c:v libvpx-vp9 -b:v 1M -c:a libopus -vf \"setpts=2.0*PTS\" "
            << boost::filesystem::absolute(static_cast<path>(visualization_path_[run-1])).append("simulation.webm");
    if (parameters_.output_video) abm::util::executeShellCommand(command.str());

    auto filepath = static_cast<path>(visualization_path_[run-1]).append("simulation.webm");
    auto rel_filepath = abm::util::getLastNFolders(filepath.string(), 4);
    site.addOutputPath(rel_filepath); //send filepath to frontend api
    site.addOutputCommand(command.str()); //send command for output video back to frontend api
}

