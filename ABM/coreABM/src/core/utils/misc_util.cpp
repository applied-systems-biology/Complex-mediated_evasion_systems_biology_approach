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

#include <dirent.h>
#include <set>
#include <sstream>
#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <bits/basic_string.h>

#include "core/utils/misc_util.h"
#include "core/basic/SphericCoordinate3D.h"
#include "core/simulation/Agent.h"
#include "core/utils/macros.h"
#include "external/json.hpp" // Include the nlohmann/json.hpp header

using json = nlohmann::json;

namespace abm::util {

    std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents) {
        std::set<std::string> sorted_content;
        std::transform(agents.begin(),
                       agents.end(),
                       std::inserter(sorted_content, sorted_content.begin()),
                       [current_time](const auto &agent) {
                           return std::to_string(agent->getId()) + "-" +
                                  std::to_string(agent->getSurface()->getAllSpheresOfThis().size()) +
                                  "-" + agent->getPosition().printCoordinates();
                       });

        std::ostringstream buff_string;
        std::copy(sorted_content.begin(), sorted_content.end(), std::ostream_iterator<std::string>(buff_string, " "));
        return std::to_string(std::hash<std::string>{}(buff_string.str()));
    }

    template<char T, typename... Args>
    std::string concatenate(Args &&... args) {
        std::ostringstream output;
        ((output << std::forward<Args>(args) << T), ...);
        return output.str();
    }

    std::vector<std::string> getFileNamesFromDirectory(const std::string &path, const std::string &fileMask) {
        std::vector<std::string> fileList;
        DIR *dirp = opendir(path.c_str());
        struct dirent *dp;
        while ((dp = readdir(dirp)) != nullptr) {
            if (static_cast<std::string>(dp->d_name).find(fileMask) != std::string::npos) {
                fileList.emplace_back(dp->d_name);
            }
        }
        std::sort(fileList.begin(),
                  fileList.end()); // apply some order on files to have equal vector on different machines
        closedir(dirp);
        return fileList;
    }

    bool folderExists(std::string folder) {
        return std::filesystem::exists(folder);
    }

    void read3DCoordinatesFromFile(std::vector<Coordinate3D> &AMpos, const std::string &inputString) {
        std::ostringstream am_dist_path;
        am_dist_path << inputString;
        std::ifstream file(am_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        while (getline(file, currentLine)) {
            std::vector<std::string> tokens;
            boost::algorithm::split(tokens, currentLine, boost::is_any_of(","));
            double x = atof(tokens.at(1).c_str());
            double y = atof(tokens.at(2).c_str());
            double z = atof(tokens.at(3).c_str());

            Coordinate3D receivedCoordinate{x, y, z};
            AMpos.emplace_back(receivedCoordinate);
        }
        file.close();
    }

    double readLambdaValueFromFile(const std::string &inputString) {
        std::ostringstream AM_dist_path;
        AM_dist_path << inputString;
        std::ifstream file(AM_dist_path.str().c_str());
        std::string currentLine;
        getline(file, currentLine);
        return std::stod(currentLine);
    }

    std::string getLastNFolders(const std::string& dir_path, const int n) {
        boost::filesystem::path path(dir_path);
        std::deque<std::string> folders;

        for (const auto& folder : path) {
            folders.push_back(folder.string());
            if (folders.size() > n) {
                folders.pop_front();
            }
        }

        std::string result;
        for (const auto& folder : folders) {
            result += "/";
            result += folder;
        }

        return result;
    }

    void writePairVectorToJsonFile(const std::vector<std::pair<std::string, std::vector<std::string>>>& vec, const std::string& filename) {
        // Create a JSON object and fill it with the elements of the vector
        json j = json::object();
        for (const auto& pair : vec) {
            j[pair.first] = pair.second;
        }

        // Open the file and write the JSON object to it
        std::ofstream file(filename);
        if (file.is_open()) {
            file << j.dump(4); // The '4' here is for pretty-printing with 4 spaces of indentation
            file.close();
        } else {
            std::cerr << "Unable to open file: " << filename << std::endl;
        }
    }

    std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv) {

        std::unordered_map<std::string, std::string> inputArgs;

        if (argc > 2) {
            if (argc % 2 == 0) {
                for (int i = 2; i < argc; i = i + 2) {
                    std::ostringstream ss, ssVal;
                    ss << argv[i];
                    std::string curArg = ss.str();
                    std::string firstCharArg = curArg.substr(0, 1);
                    if (firstCharArg == "-") {
                        ssVal << argv[i + 1];
                        std::string curVal = ssVal.str();
                        std::string withoutFirstCharArg = curArg.substr(1, -1);
                        inputArgs[withoutFirstCharArg] = curVal;
                    } else {
                        ERROR_STDERR(
                                R"(usage hint: "./hABM <config.json> -id1 value1 -id2 value2 ... ")");
                        ERROR_STDERR("now stopping execution");
                        exit(1);
                    }
                }
            } else {
                ERROR_STDERR(
                        R"(usage hint: "./hABM <config.json> -id1 value1 -id2 value2 ... ")");
                ERROR_STDERR("now stopping execution");
                exit(1);
            }
        }

        return inputArgs;
    }

    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para) {
        std::vector<std::string> keys;
        std::vector<std::vector<std::string>> values;
        std::vector<int> its;
        std::vector<std::vector<std::string>> tuples;
        for (const auto &x: para) {
            keys.emplace_back(x.first);
            values.emplace_back(x.second);
            its.emplace_back(0);
        }

        bool keep_running = true;
        while (keep_running) {
            std::vector<std::string> valtuple;
            for (int i = 0; i < keys.size(); ++i) {
                valtuple.emplace_back(values[i][its[i]]);
            }
            tuples.emplace_back(valtuple);

            for (int j = 0; j < keys.size();) {
                if (its[j] < values[j].size() - 1) {
                    its[j] += 1;
                    j = keys.size();
                } else {
                    if (j == keys.size() - 1) {
                        keep_running = false;
                    }
                    its[j] = 0;
                    ++j;
                }
            }

        }
        return {keys, tuples};
    }

    bool approxEqual(double d1, double d2, double epsilon) {
        return std::abs(d1 - d2) < epsilon;
    }

    Coordinate3D toCartesianCoordinates(const SphericCoordinate3D &coordinate) {
        double x = coordinate.r * sin(coordinate.theta) * cos(coordinate.phi);
        double y = coordinate.r * sin(coordinate.theta) * sin(coordinate.phi);
        double z = coordinate.r * cos(coordinate.theta);

        return {x, y, z};
    }

    SphericCoordinate3D toSphericCoordinates(const Coordinate3D &coordinate) {
        double r, theta, phi;
        r = sqrt(coordinate.x * coordinate.x + coordinate.y * coordinate.y + coordinate.z * coordinate.z);
        if (r == 0) {
            return {0, 0, 0};
        }
        theta = acos(coordinate.z / r);
        if (coordinate.x > 0) {
            phi = atan(coordinate.y / coordinate.x);
        } else {
            if (coordinate.x == 0.0) {
                if (coordinate.y >= 0) {
                    phi = M_PI / 2;
                } else {
                    phi = -M_PI / 2;
                }
            } else {
                if (coordinate.y >= 0.0) {
                    phi = atan(coordinate.y / coordinate.x) + M_PI;
                } else {
                    phi = atan(coordinate.y / coordinate.x) - M_PI;
                }
            }
        }
        return {r, theta, phi};
    }

    Coordinate3D rotateAboutThetaAndPhi(const Coordinate3D &coordinate, double theta, double phi) {
        SphericCoordinate3D pos_sph_old = toSphericCoordinates(coordinate);
        return toCartesianCoordinates({pos_sph_old.r, pos_sph_old.theta + theta, pos_sph_old.phi + phi});
    }

    double dotProduct(const Coordinate3D &p, const Coordinate3D &q) {
        return p.x * q.x + p.y * q.y + p.z * q.z;
    }

    double angleBetweenVectors(const Coordinate3D &p, const Coordinate3D &q) {
        double dotp = dotProduct(p, q);
        double magp = p.getMagnitude();
        double magq = q.getMagnitude();
        return acos(dotp / (magp * magq));
    }

    double angleThreePoints(const Coordinate3D &p, const Coordinate3D &q, const Coordinate3D &u) {
        Coordinate3D vec1 = p - q;
        Coordinate3D vec2 = u - q;
        return angleBetweenVectors(vec1, vec2);
    }

    Coordinate3D genPerpendicularVector(const Coordinate3D &n) {
        double x = 1.0, y = 1.0, z = 1.0;
        if (n.x != 0) {
            x = (-n.y * y - n.z * z) / n.x;
        } else if (n.y != 0) {
            y = (-n.x * x - n.z * z) / n.y;
        } else {
            z = (-n.x * x - n.y * y) / n.z;
        }
        Coordinate3D res = {x, y, z};
        return res;
    }

    Coordinate3D genPerpendicularVector(const Coordinate3D &n, const Coordinate3D &m) {
        double x = 1.0, y = 1.0, z = 1.0;
        if (m.y * n.x - n.y * m.x != 0) {
            if (n.x != 0) {
                y = (m.x * (n.z * z) + n.x * (-m.z * z)) / (m.y * n.x - n.y * m.x);
                x = (-n.y * y - n.z * z) / n.x;
            } else if (m.x != 0) {
                y = (m.x * (n.z * z) + n.x * (-m.z * z)) / (m.y * n.x - n.y * m.x);
                x = (-m.y * y - m.z * z) / m.x;
            } else {
                ERROR_STDERR("Error in vector calculation; denominator is zero");
            }
        } else {
            ERROR_STDERR("Error in vector calculation; denominator is zero");
        }
        Coordinate3D res = {x, y, z};

        if (!(abm::util::approxEqual(dotProduct(n, res), 0.0) || abm::util::approxEqual(dotProduct(m, res), 0.0))) {
            ERROR_STDERR("Error in vector calculation; vectors are not perpendicular");
        }
        return res;
    }

    std::pair<double, double> getRandom2Direction(Randomizer* randomizer, double mean, double std) {
        std::pair<double, double> direction;
        double angle = 0;
        if (mean == 0 && std == 0) {
            angle = randomizer->generateDouble(0, M_PI);
        } else {
            angle = randomizer->generateNormalDistributedValue(mean, std);
        }

        if (abm::util::approxEqual(angle, M_PI)) {
            direction.first = 0.0;
            direction.second = (randomizer->generateInt(0, 1) - 0.5);
        } else {
            if (angle > M_PI / 2) {
                direction.first = -1.0;
            } else {
                direction.first = 1.0;
            }
            direction.second = (randomizer->generateInt(0, 1) - 0.5) * 2 * abs(tan(angle));
        }

        return direction;
    }

    std::vector<unsigned int> generateRandomPermutation(Randomizer *randomizer, unsigned int size) {
        std::vector<unsigned int> permutationVector;
        if (size > 0) {
            //initialize the vector
            for (unsigned int k = 0; k < size; k++) {
                permutationVector.push_back(k);
            }
            for (unsigned int k = 0; k < size - 1; k++) {
                swap(&permutationVector, k, randomizer->generateInt(k, size - 1));
            }
        }
        return permutationVector;
    }

    void swap(std::vector<unsigned int> *vectorToSwap, unsigned int i, unsigned int j) {
        unsigned int temp;
        temp = (*vectorToSwap)[i];
        (*vectorToSwap)[i] = (*vectorToSwap)[j];
        (*vectorToSwap)[j] = temp;
    }

    bool isSubstring(const std::string& str1, const std::string& str2) {
        if (str1.size() > str2.size()) {
            if (str1.find(str2) != std::string::npos) {
                return true;
            } else {
                return false;
            }
        } else {
            if (str2.find(str1) != std::string::npos) {
                return true;
            } else {
                return false;
            }
        }
    }

    double nChoosek(unsigned long n, unsigned long k) {
        if (k * 2 > n) { k = n - k; }

        if (k > n) return 0;
        if (k == 0) return 1;

        double result = 1.0;
        for (unsigned long i = 1; i <= k; ++i) {
            result *= (n - i + 1);
            result /= i;
        }
        return result;
    }

    double bernoulliProbability(unsigned long n, unsigned long k, double p) {
        return nChoosek(n, k) * pow(p, k) * pow(1 - p, n - k);
    }

    bool compareTriple(const std::vector<unsigned int>& a, const std::vector<unsigned int>& b) {
        return std::tie(a[0], a[1], a[2]) < std::tie(b[0], b[1], b[2]);
    }
}