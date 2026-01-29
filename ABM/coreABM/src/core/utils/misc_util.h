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

#ifndef CORE_UTILS_MISC_UTIL_H
#define CORE_UTILS_MISC_UTIL_H

#include <vector>
#include <string>
#include <memory>
#include <unordered_map>

#include "core/basic/Coordinate3D.h"
#include "core/basic/SphericCoordinate3D.h"

class Agent;
class Randomizer;
namespace abm::util {
    template<typename... Ts>
    struct overloaded : Ts ... {
        using Ts::operator()...;
    };
    template<typename... Ts> overloaded(Ts...) -> overloaded<Ts...>;
    template<char T = ';', typename... Args>
    std::string concatenate(Args &&... args);

    /*!
     * Generates hash from positions of outputs to compare simulation results
     * @param current_time Double for current time
     * @param agents vector of Agent that contains agents of current simulations
     * @return String that contains calculated hash value
     */
    std::string generateHashFromAgents(const double &current_time, const std::vector<std::shared_ptr<Agent>> &agents);

    /*!
     * Get file names from a directory
     * @param path String that contains directory
     * @param fileMask String that contains a file mask to only get files of one type (i.e. *.csv)
     * @return vector of String that contains all file names
     */
    std::vector<std::string> getFileNamesFromDirectory(const std::string &path, const std::string &fileMask = "");

    /*!
     * Generates cartesian product of n input sets
     * @param para Unordered map of Strings of parameters and corresponding values
     * @return pair of vector of String that contains all combinations of given parameters values
     */
    std::pair<std::vector<std::string>, std::vector<std::vector<std::string>>>
    calculateCartesianProd(const std::unordered_map<std::string, std::vector<std::string>> &para);

    /*!
     * Read coordinates from file
     * @param AMpos vector of Coordinate3D that contains written positions from file
     * @param inputString String that contains path to file
     */
    void read3DCoordinatesFromFile(std::vector<Coordinate3D> &AMpos, const std::string &inputString);

    /*!
     * read Lambda value from file for input of AM
     * @param input_string String that contains path to file
     * @return Double that contains Lambda value for input of AM
     */
    double readLambdaValueFromFile(const std::string &input_string);

    /*!
     * Checks if a folder exists
     * @return String that contains path to folder
     */
    bool folderExists(std::string);

    // get last n folders of string path
    std::string getLastNFolders(const std::string& dir_path, const int n);

    // writes vector of pairs to json file
    void writePairVectorToJsonFile(const std::vector<std::pair<std::string, std::vector<std::string>>>& vec, const std::string& filename);

    /*!
     * Checks if two values are approximately equal up to a certain accuracy
     * @param d1 Double that contains first number
     * @param d2 Double that contains second number
     * @param epsilon Double that contains accuracy
     * @return Bool if values are equal or not
     */
    bool approxEqual(double d1, double d2, double epsilon = 1e-8);

    /// Calculates dotProduct
    double dotProduct(const Coordinate3D &p, const Coordinate3D &q);

    /// Calculates angle between vectors
    double angleBetweenVectors(const Coordinate3D &p, const Coordinate3D &q);

    /// Calculates angles between three points
    double angleThreePoints(const Coordinate3D &p, const Coordinate3D &q, const Coordinate3D &u);

    /// Generates a perpendicular vector to one input vector
    Coordinate3D genPerpendicularVector(const Coordinate3D &n);

    /// Generates a perpendicular vector to two input vectors
    Coordinate3D genPerpendicularVector(const Coordinate3D &n, const Coordinate3D &m);

    /// Generates a random direction in 2D
    std::pair<double, double> getRandom2Direction(Randomizer* randomizer, double mean=0.0, double std=0.0);

    /// Transforms SphericCoordinate3D to Coordinate3D
    Coordinate3D toCartesianCoordinates(const SphericCoordinate3D &coordinate);

    /// Transforms Coordinate3D to SphericCoordinate3D
    SphericCoordinate3D toSphericCoordinates(const Coordinate3D &coordinate);

    /// Rotates a Coordinate3D about an angle theta and phi
    Coordinate3D rotateAboutThetaAndPhi(const Coordinate3D &coordinate, double theta, double phi);

    /// Handles cmd inputs
    std::unordered_map<std::string, std::string> handleCmdInputs(int argc, char **argv);

    std::vector<unsigned int> generateRandomPermutation(Randomizer *randomizer, unsigned int size);

    void swap(std::vector<unsigned int> *vectorToSwap, unsigned int i, unsigned int j);

    bool isSubstring(const std::string& str1, const std::string& str2);

    double nChoosek(unsigned long n, unsigned long k);

    double bernoulliProbability(unsigned long n, unsigned long k, double p);

    bool compareTriple(const std::vector<unsigned int>& a, const std::vector<unsigned int>& b);
}
#endif //CORE_UTILS_MISC_UTIL_H
