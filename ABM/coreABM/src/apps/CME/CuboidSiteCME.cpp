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

#include "CuboidSiteCME.h"

#include "InSituMeasurementsCME.h"
#include "apps/CME/cells/Complex.h"
#include "apps/CME/cells/Drug.h"
#include "apps/CME/cells/Pathogenic_cell.h"
#include "apps/CME/factories/CellFactoryCME.h"
#include "apps/CME/factories/InteractionFactoryCME.h"
#include "apps/CME/factories/InteractionStateFactoryCME.h"
#include "core/simulation/AgentManager.h"
#include "core/simulation/factories/CellFactory.h"
#include "core/simulation/neighbourhood/BalloonListNHLocator.h"
#include "core/simulation/neighbourhood/NeighbourhoodLocator.h"
#include "io_utilsCME.h"
#include <utility>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>

#include "factories/RateFactoryCME.h"

namespace fs = std::filesystem;


CuboidSiteCME::CuboidSiteCME(Randomizer* random_generator,
                                   std::shared_ptr<InSituMeasurements> measurements,
                                   std::string config_path,
                                   std::unordered_map<std::string, std::string> cmd_input_args,
                                   std::string output_path)
                                   : Site(random_generator, std::move(measurements), config_path,  cmd_input_args, output_path) {

    //SYSTEM_STDOUT("Initialize CuboidSiteCME ...");
    measurements_->setSite(this);
    const auto simulator_config = static_cast<boost::filesystem::path>(config_path).append("simulator-config.json");
    parameters_ = abm::utilCME::getSimulationParameters(simulator_config.string());

    // Read in input parameters
    const auto input_config = static_cast<boost::filesystem::path>(config_path).append("input-config.json");
    input_parameters_ = abm::utilCME::getInputParameters(input_config.string());

    // Read in cmd input and screening parameters
    parameters_.cmd_input_args = cmd_input_args;
    if (!parameters_.cmd_input_args.empty()) {
        handleCmdInputArgs(cmd_input_args);
    }

    // Initialize factories
    interaction_factory_ = std::make_unique<InteractionFactoryCME>(parameters_.interaction_parameters, parameters_.use_interactions);
    rate_factory_ = std::make_unique<RateFactoryCME>(input_parameters_.rates, this->getTimeStepping());
    cell_factory_ = std::make_unique<CellFactoryCME>(parameters_.site_parameters);
    cell_state_factory_ = std::make_unique<CellStateFactory>(parameters_.site_parameters, rate_factory_.get());
    interaction_state_factory_ = std::make_unique<InteractionStateFactoryCME>(parameters_.interaction_parameters, rate_factory_.get());
    agent_manager_ = std::make_unique<AgentManager>(parameters_, this);

    dimensions = parameters_.dimensions;
    //setStoppingCondition(parameters_.stopping_criteria);

    auto* cuboid_parameters = static_cast<abm::util::SimulationParameters::CuboidSiteParameters*>(parameters_.site_parameters.get());
    identifier_ = cuboid_parameters->identifier;
    upper_bound_ = cuboid_parameters->upper_bound;
    lower_bound_ = cuboid_parameters->lower_bound;
    molecules_grid_size_ = cuboid_parameters->molecules_grid_size;
    boundary_type_ = cuboid_parameters->boundary_condition;
    setBoundaryCondition();
    neighbourhood_locator_ = std::make_unique<BalloonListNHLocator>(cuboid_parameters->agent_manager_parameters.agents, lower_bound_, upper_bound_, this);
    initializeAgents(cuboid_parameters->agent_manager_parameters, 0, parameters_.time_stepping);
    agent_manager_->setInitFungalQuantity();

    auto* cuboid_CME_parameters = static_cast<abm::utilCME::CuboidSiteCMEParameters*>(parameters_.site_parameters.get());

    flow_ = cuboid_CME_parameters->flow;

    for (const auto& name : getAgentManager()->getAllAgentTypes()) {
        const auto nb = getAgentManager()->getAgentQuantity(name);
        prev_populations[name] = nb;
    }

    stopping_criteria_ = cuboid_CME_parameters->stopping_criteria;
    stop_sim = false;
}

bool CuboidSiteCME::containsPosition(Coordinate3D position) {
    bool contains = (position.x >= lower_bound_.x) && (position.x <= upper_bound_.x) &&
                    (position.y >= lower_bound_.y) && (position.y <= upper_bound_.y) &&
                    (position.z >= lower_bound_.z) && (position.z <= upper_bound_.z);

    return contains;
}

Coordinate3D CuboidSiteCME::getRandomPosition(double cell_radius) {
    double x, y, z;
    double dist;
    bool insideSphere;

    do {
        x = random_generator_->generateDouble(upper_bound_.x, lower_bound_.x);
        y = random_generator_->generateDouble(upper_bound_.y, lower_bound_.y);
        z = random_generator_->generateDouble(upper_bound_.z, lower_bound_.z);

        auto pos = Coordinate3D{x, y, z};
        dist = pos.calculateEuclidianDistance(Coordinate3D{0, 0, 0});
        insideSphere = (dist < cell_radius);
    } while (insideSphere);

    return Coordinate3D{x, y, z};
}

void CuboidSiteCME::handleBoundaryCross(Agent* agent, Coordinate3D* move_vec, double current_time) {

    boundary_condition_->handleBoundaryCross(agent, move_vec, current_time);
}

Coordinate3D CuboidSiteCME::getLowerLimits() {
    return lower_bound_;
}

Coordinate3D CuboidSiteCME::getUpperLimits() {
    return upper_bound_;
}

Coordinate3D CuboidSiteCME::getRandomBoundaryPoint() {
    double x, y, z;
    double ylength, xlength;
    double xProb;

    double xyArea, xzArea, yzArea, totalArea;
    double xyProb, xzProb;

    double rand;
    unsigned int use_lower_bound_;

    xyArea = fabs((upper_bound_.y - lower_bound_.y) * (upper_bound_.x - lower_bound_.x));
    xzArea = fabs((upper_bound_.x - lower_bound_.x) * (upper_bound_.z - lower_bound_.z));
    yzArea = fabs((upper_bound_.y - lower_bound_.y) * (upper_bound_.z - lower_bound_.z));

    totalArea = xyArea + xzArea + yzArea;

    xyProb = xyArea / totalArea;
    xzProb = xzArea / totalArea;

    rand = random_generator_->generateDouble();
    use_lower_bound_ = random_generator_->generateInt(1);
    if (rand < xyProb) {
        //place point in xyPlane
        if (use_lower_bound_) {
            x = random_generator_->generateDouble(lower_bound_.x, upper_bound_.x);
            y = random_generator_->generateDouble(lower_bound_.y, upper_bound_.y);
            z = lower_bound_.z;
        } else {
            x = random_generator_->generateDouble(lower_bound_.x, upper_bound_.x);
            y = random_generator_->generateDouble(lower_bound_.y, upper_bound_.y);
            z = upper_bound_.z;
        }
    } else {
        if (rand < xyProb + xzProb) {
            //place point in xzPlane
            if (use_lower_bound_) {
                x = random_generator_->generateDouble(lower_bound_.x, upper_bound_.x);
                y = lower_bound_.y;
                z = random_generator_->generateDouble(lower_bound_.z, upper_bound_.z);
            } else {
                x = random_generator_->generateDouble(lower_bound_.x, upper_bound_.x);
                y = upper_bound_.y;
                z = random_generator_->generateDouble(lower_bound_.z, upper_bound_.z);
            }
        } else {
            //place point in yzPlane
            if (use_lower_bound_) {
                x = lower_bound_.x;
                y = random_generator_->generateDouble(lower_bound_.y, upper_bound_.y);
                z = random_generator_->generateDouble(lower_bound_.z, upper_bound_.z);
            } else {
                x = upper_bound_.x;
                y = random_generator_->generateDouble(lower_bound_.y, upper_bound_.y);
                z = random_generator_->generateDouble(lower_bound_.z, upper_bound_.z);
            }
        }
    }


    if (dimensions == 2) {
        z = (upper_bound_.z + lower_bound_.z)/2;
    }
    Coordinate3D nextPos = Coordinate3D();
    Coordinate3D boundaryPoint = Coordinate3D{x, y, z};

    do {
        boundary_input_vector_ = generateRandomDirectionVector(boundaryPoint, 1);
        nextPos = boundaryPoint;
        nextPos += boundary_input_vector_;
    } while (!containsPosition(nextPos));

    return Coordinate3D{x, y, z};
}

Coordinate3D CuboidSiteCME::generateRandomDirectionVector(Coordinate3D position, double length) {

    //get a sin() sampled value in [0,PI] for theta via a uniform value of u
    if (dimensions == 2) {
        double phi = random_generator_->generateDouble(M_PI * 2.0);
        double x = length * cos(phi);
        double y = length * sin(phi);
        double z = 0;
        return Coordinate3D{x, y, z};
    } else if (dimensions == 3) {
        double u = random_generator_->generateDouble(); //sampler->sample();

        double phi = random_generator_->generateDouble(M_PI * 2.0);
        double r = length;

        double subst = 2 * r * sqrt(u * (1 - u));
        double x = subst * cos(phi);
        double y = subst * sin(phi);
        double z = r * (1 - 2 * u);
        return Coordinate3D{x, y, z};
    } else {
        return Coordinate3D{};
    }
}

Coordinate3D CuboidSiteCME::generatePersistentDirectionVector(Coordinate3D position,
                                                           double length,
                                                           Coordinate3D prevVector,
                                                           double previousAlpha) {
    if (prevVector.getMagnitude() == 0){
        prevVector *= 0;
    }
    else {
        prevVector *= length/prevVector.getMagnitude();
    }
    return prevVector;
}

Coordinate3D CuboidSiteCME::generateBackShiftOnContacting(SphereRepresentation* activeSphere, SphereRepresentation* passiveSphere, double mustOverhead) {

    Coordinate3D normalVecBetwCells = activeSphere->getEffectiveConnection(passiveSphere);
    double overlap = passiveSphere->getRadius() + activeSphere->getRadius() - normalVecBetwCells.getMagnitude();

    Coordinate3D backShift = normalVecBetwCells;
    if (overlap > 0 && backShift.getMagnitude() !=0 ) {
        backShift *= -(overlap - mustOverhead)/backShift.getMagnitude();
    } else {
        backShift *=0;
    }

    return backShift;
}

Coordinate3D CuboidSiteCME::getCenterPosition() {

    double   x = (upper_bound_.x + lower_bound_.x) / 2.0;
    double   y = (upper_bound_.y + lower_bound_.y) / 2.0;
    double   z = (upper_bound_.z + lower_bound_.z) / 2.0;

    return Coordinate3D{x, y, z};
}

void CuboidSiteCME::initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                                  double current_time,
                                  double time_delta) {
    auto count_cell = 0;
    auto cell_radius = 0.0;
    for (const auto &agent_parameters: parameters.agents) {
        if (agent_parameters->type == "Pathogenic_cell") {
            cell_radius = agent_parameters->morphology_parameters.radius;
            break;
        }
    }
    for (const auto &agent_parameters: parameters.agents) {
        const auto init_distribution = agent_parameters->initial_distribution;
        const auto name = agent_parameters->type;
        const auto number_of_agents = agent_parameters->number;
        input_rates_[name] = agent_parameters->input_lambda;
        double z_level = (upper_bound_.z + lower_bound_.z)/2.0;
        for (auto i = 0; i < number_of_agents; ++i) {
            auto radius = agent_parameters->morphology_parameters.radius;
            Coordinate3D initial_position;
            switch (init_distribution) {
                case 0:
                    initial_position = getRandomPosition(cell_radius);
                    break;
                case 1:
                    initial_position = getCenterPosition();
                    break;
                case 2:
                    initial_position = random_point_on_sphere(cell_radius, Coordinate3D{0,0,0});
                    break;
                case 3:
                    do {
                        initial_position = getZEqualsZeroPosition(agent_parameters->morphology_parameters.radius);
                    } while (check_collision(initial_position));
                    break;
                case 4:
                    ++count_cell;
                    initial_position = two_cells_layout(radius, count_cell, 3);
                    break;
                default:
                    initial_position = getRandomPosition(cell_radius);
                    break;
            }
            Coordinate3D initial_vector = Coordinate3D();
            if (dimensions == 2) { initial_position.z = z_level;}
            auto cell_complex = "";
            if (agent_parameters->type == "Complex") {
                cell_complex = "AMP";
            }
            auto agent = agent_manager_->emplace_back(cell_factory_->createCell(name,
                                                                              std::make_unique<Coordinate3D>(initial_position),
                                                                              agent_manager_->generateNewID(),
                                                                              this,
                                                                              time_delta,
                                                                              current_time,
                                                                              cell_complex));
            if (initial_vector.getMagnitude() != 0) {
                agent->getMovement()->setPreviousMove(&initial_vector);
            }
            if (agent->getTypeName() == "FungalCell") {
                agent_manager_->addFungalCellToList(agent.get());
            }
        }
    }
}

Coordinate3D CuboidSiteCME::random_point_on_sphere(double radius, Coordinate3D pos){
    double x = random_generator_->generateDouble(-1.0, 1.0);
    double y = random_generator_->generateDouble(-1.0, 1.0);
    double z = random_generator_->generateDouble(-1.0, 1.0);

    double length = std::sqrt(x*x + y*y + z*z);
    auto res = Coordinate3D{x, y, z};
    res *= radius/length;
    res += pos;

    return res;
}

Coordinate3D CuboidSiteCME::biased_random_point_on_sphere(double radius, Coordinate3D pos_of_cell, Coordinate3D pos_of_uptake){
    // How do we sample a biased random pos on radius based on a given pos?
    double x = random_generator_->generateDouble(-1.0, 1.0);
    double y = random_generator_->generateDouble(-1.0, 1.0);
    double z = random_generator_->generateDouble(-1.0, 1.0);

    double length = std::sqrt(x*x + y*y + z*z);
    auto res = Coordinate3D{x, y, z};
    res *= radius/length;
    res += pos_of_cell;

    // Bias towards pos_of_uptake
    Coordinate3D bias = pos_of_uptake - pos_of_cell;
    double bias_factor = 0.9; // Adjust bias factor as needed
    res += bias * bias_factor;

    double res_length = std::sqrt(res.x * res.x + res.y * res.y + res.z * res.z);
    res *= radius / res_length;

    //DEBUG_STDOUT("Dist new position = " << res.calculateEuclidianDistance(pos_of_cell));
    return res;
}

void CuboidSiteCME::secretion_of_agents_at_cell_surface_randomly(double timestep, double current_time, double rate, const std::string& agent_type, Coordinate3D pos, double radius) {
    auto tmp = rate * timestep;// + secretion_record_[agent_type];
    auto nb_mol = floor(tmp);
    auto proba = tmp - nb_mol;
    double const p = this->random_generator_->generateDouble();
    for (int i = 0; i<nb_mol; ++i){
        agent_manager_->emplace_back(cell_factory_->createCell(agent_type, std::make_unique<Coordinate3D>(random_point_on_sphere(radius, pos)),
                                                               agent_manager_->generateNewID(), this, timestep, current_time, ""));
    }
    if (p < proba){
        agent_manager_->emplace_back(cell_factory_->createCell(agent_type, std::make_unique<Coordinate3D>(random_point_on_sphere(radius, pos)),
                                                               agent_manager_->generateNewID(), this, timestep, current_time, ""));
    }
}

void CuboidSiteCME::secretion_of_agents_at_cell_surface(double timestep, double current_time, double rate, const std::string& agent_type, Coordinate3D pos_of_cell,double radius, std::vector<Coordinate3D> pos_of_uptake) {
    auto tmp = rate * timestep;// + secretion_record_[agent_type];
    auto nb_mol = floor(tmp);
    auto proba = tmp - nb_mol;
    double const p = this->random_generator_->generateDouble();
    for (int i = 0; i<nb_mol; ++i){
        agent_manager_->emplace_back(cell_factory_->createCell(agent_type, std::make_unique<Coordinate3D>(biased_random_point_on_sphere(radius, pos_of_cell, pos_of_uptake[i])),
                                                               agent_manager_->generateNewID(), this, timestep, current_time, ""));
    }
    if (p < proba){
        auto pos = biased_random_point_on_sphere(radius, pos_of_cell, pos_of_uptake.back());
        agent_manager_->emplace_back(cell_factory_->createCell(agent_type, std::make_unique<Coordinate3D>(pos),
                                                               agent_manager_->generateNewID(), this, timestep, current_time, ""));
    }
}


void CuboidSiteCME::doAgentDynamics(Randomizer *random_generator, SimulationTime &time) {
    std::static_pointer_cast<InSituMeasurementsCME>(measurements_)->observeMeasurements(time);

    const auto current_time = time.getCurrentTime();
    const auto dt = time.getCurrentDeltaT();

    const auto &all_agents = agent_manager_->getAllAgents();
    if (!all_agents.empty()) {
        // Loop over all agents (random order)
        const auto current_order = abm::util::generateRandomPermutation(random_generator, all_agents.size());
        for (auto agent_idx = current_order.begin(); agent_idx < current_order.end(); ++agent_idx) {
            auto curr_agent = all_agents[*agent_idx];
            if (nullptr != curr_agent) {
                // Do all actions for one timestep for each agent (-> Cell.cpp)
                curr_agent->doAllActionsForTimestep(dt, current_time);
            }
        }
        // Clean up agents
        agent_manager_->cleanUpAgents(current_time);
        //agent_manager_->inputOfAgents(current_time, random_generator);
    }
    if (!stop_sim) {
        stop_sim = check_steady_state(time.getCurrentDeltaT());
    }
}

Coordinate3D CuboidSiteCME::getZEqualsZeroPosition(double radius) {

    double   x = random_generator_->generateDouble(lower_bound_.x + radius, upper_bound_.x - radius);
    double   y = random_generator_->generateDouble(lower_bound_.y + radius, upper_bound_.y - radius);

    return Coordinate3D{x, y, 0.0};
}

bool CuboidSiteCME::check_collision(Coordinate3D pos_to_check) {
    bool check = false;
    for (const auto &ag: agent_manager_->getAllAgents()){
        if (ag->getPosition().calculateEuclidianDistance(pos_to_check) <= ag->getMorphology()->getBasicSphereOfThis()->getRadius()*2) {
            check = true;
            break;
        }
    }
    return check;
}

void CuboidSiteCME::handleCmdInputArgs(std::unordered_map<std::string, std::string>cmd_input_args) {
    if (cmd_input_args.size() > 0) {
        for (const auto &[key, value]: cmd_input_args) {
            for (auto &input_rate: input_parameters_.rates) {
                if (abm::util::isSubstring(key,input_rate->key)) {
                    if (abm::util::isSubstring(key, "alpha")) {
                        auto conc_rate = dynamic_cast<abm::utilCME::ConcentrationDependentRateParameters*>(input_rate.get());
                        conc_rate->alpha = std::stod(value);
                        DEBUG_STDOUT("conc_rate->alpha = " << conc_rate->alpha);
                    }
                    else {
                        input_rate->rate = std::stod(value);
                        DEBUG_STDOUT("input_rate->rate = " << input_rate->rate);
                    }
                }
            }
            if (abm::util::isSubstring(key, "secretion")){
                int char_pos = key.find("_");
                auto mol = key.substr(char_pos+1);
                auto& ag_param = parameters_.site_parameters->agent_manager_parameters;
                for (auto& ag: ag_param.agents){
                    if (ag->type == "Pathogenic_cell") {
                        auto ag_sec = std::static_pointer_cast<abm::utilCME::AgentsParametersCME>(ag);
                        for (auto& [name, data] : ag_sec->secretion_) {
                            if (name == mol) {
                                data.rate = std::stod(value);
                            }
                        }
                    }
                }
            }
            if (abm::util::isSubstring(key, "lag")) {
                int char_pos = key.find("_");
                auto mol = key.substr(char_pos + 1);
                auto& ag_param = parameters_.site_parameters->agent_manager_parameters;
                for (auto& ag : ag_param.agents) {
                    if (ag->type == "Pathogenic_cell") {
                        auto ag_sec = std::static_pointer_cast<abm::utilCME::AgentsParametersCME>(ag);
                        for (auto& [name, data] : ag_sec->secretion_) {
                            if (name == mol) {
                                data.lag = std::stod(value);
                            }
                        }
                    }
                }
            }
            if (key == "nb_receptor_blocked") {
                // Define relative file path using filesystem::path
                fs::path file_path = fs::path("configurations") / "basicConfigs" / "basicConfigCME" / "input" / "receptor_radius_diffusion_downscaled_MW_based.csv";
                auto& ag_param = parameters_.site_parameters->agent_manager_parameters;
                for (auto& ag : ag_param.agents) {
                    if (ag->type == "Drug") {
                        auto ag_drug = std::static_pointer_cast<abm::utilCME::Drug>(ag);
                        ag_drug->nb_recep_blocked = std::stod(value);
                        SYSTEM_STDOUT("nb_recep_blocked = " << ag_drug->nb_recep_blocked);
                        // Open file inside loop (ensures fresh read for each agent)
                        std::ifstream file(file_path);
                        if (!file.is_open()) {
                            throw std::runtime_error("Could not open file: " + file_path.string());
                        }

                        std::string line, header;
                        std::getline(file, header); // Skip the header line

                        bool match_found = false;
                        while (std::getline(file, line)) {
                            std::stringstream ss(line);
                            std::string cell;
                            std::vector<std::string> row;

                            // Read CSV row into vector
                            while (std::getline(ss, cell, ',')) {
                                row.push_back(cell);
                            }

                            int nb_receptors = std::stoi(row[0]);
                            if (ag_drug->nb_recep_blocked == nb_receptors) {
                                // Assign values
                                ag_drug->morphology_parameters.radius = std::stod(row[1]);
                                //ag_drug->movement_parameters.diffusion_coefficient = std::stod(row[2]);

                                SYSTEM_STDOUT("radius = " << ag_drug->morphology_parameters.radius);
                                //SYSTEM_STDOUT("D = " << ag_drug->movement_parameters.diffusion_coefficient);
                                match_found = true;
                                break; // Stop searching once a match is found
                            }
                        }

                        file.close(); // Close file after processing

                        if (!match_found) {
                            throw std::runtime_error("nb_receptor_blocked value not found in CSV: " + value);
                        }
                    }
                }
            }
            if (key == "drug_conc") {
                auto& ag_param = parameters_.site_parameters->agent_manager_parameters;
                for (auto& ag : ag_param.agents) {
                    if (ag->type == "Drug") {
                        auto ag_drug = std::static_pointer_cast<abm::utilCME::Drug>(ag);
                        ag_drug->number = std::stod(value);
                    }
                }
            }
            if (key == "drug_D") {
                auto& ag_param = parameters_.site_parameters->agent_manager_parameters;
                for (auto& ag : ag_param.agents) {
                    if (ag->type == "Drug") {
                        auto ag_drug = std::static_pointer_cast<abm::utilCME::Drug>(ag);
                        ag_drug->movement_parameters.diffusion_coefficient = std::stod(value);
                        SYSTEM_STDOUT("D = " << ag_drug->movement_parameters.diffusion_coefficient);
                    }
                }
            }
        }
    }
}

void CuboidSiteCME::do_site_dynamics(double timestep, double current_time) {
    for (const auto &[name, rate]: flow_) {
        double const proba = rate * timestep;
        double const p = this->random_generator_->generateDouble();
        if (p < proba){
            agent_manager_->emplace_back(cell_factory_->createCell(name, std::make_unique<Coordinate3D>(getRandomBoundaryPoint()),
                                                                   agent_manager_->generateNewID(), this, timestep, current_time, ""));
        }
    }
}
Coordinate3D CuboidSiteCME::two_cells_layout(double radius, int cell, double pos_dist) {
    double x;
    double l = upper_bound_.x - lower_bound_.x;
    double m = std::min(pos_dist, l-4*radius);
    if (m< 0){
        m = (l-4*radius)/3;
    }
    auto scale = (l-4*radius-m)/2;
    if (cell == 1){
        x = lower_bound_.x + radius +scale;
    }
    else{
        x = upper_bound_.x - scale - radius;
    }
    return Coordinate3D{x, 0.0, 0.0};
}
bool CuboidSiteCME::check_steady_state(double dt) {
    if (stopping_criteria_.molecule.empty()){
        return false;
    }
    auto name = stopping_criteria_.molecule;
    const auto current_nb = getAgentManager()->getAgentQuantity(name);
    auto diff = std::abs(current_nb - prev_populations[name]);
    if (diff < stopping_criteria_.threshold){
        stopping_criteria_.time -= dt;
    }
    else{
        prev_populations[name] = current_nb;
        stopping_criteria_.time = stopping_criteria_.init_time;
        //std::cout << "TIME RESET" << std::endl;
    }
    //std::cout << "TIME BEFORE STEADY STATE REACHED" << stopping_criteria.time << std::endl;
    return (stopping_criteria_.time <= 0);
}
