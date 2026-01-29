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

#include "CuboidSite.h"

#include <utility>
#include "core/simulation/neighbourhood/BalloonListNHLocator.h"
#include "core/simulation/neighbourhood/NeighbourhoodLocator.h"
#include "core/simulation/factories/CellFactory.h"
#include "core/simulation/AgentManager.h"

CuboidSite::CuboidSite(Randomizer* random_generator,
                       std::shared_ptr<InSituMeasurements> measurements,
                       std::string config_path,
                       std::unordered_map<std::string, std::string> cmd_input_args,
                       std::string output_path)
                       : Site(random_generator, std::move(measurements), config_path,  cmd_input_args, output_path) {

    SYSTEM_STDOUT("Initialize normal CuboidSite ...");

    // Read in simulator parameters
    auto sid = cmd_input_args.find("sid") != cmd_input_args.end() ? cmd_input_args.find("sid")->second : "";
    const auto simulator_config = static_cast<boost::filesystem::path>(config_path).append("simulator-config.json");
    parameters_ = abm::util::getSimulationParameters(simulator_config.string());

    // Read in input parameters
    const auto input_config = static_cast<boost::filesystem::path>(config_path).append("input-config.json");
    input_parameters_ = abm::util::getInputParameters(input_config.string());

    // Read in cmd inputs and screening parameters
    parameters_.cmd_input_args = cmd_input_args;
    handleCmdInputArgs(cmd_input_args);

    // Read parameters from frontend if frontend-api.json exists
    receiveFrontendParameter(parameters_, input_parameters_, config_path, output_path, sid);

    // Initialize factories
    interaction_factory_ = std::make_unique<InteractionFactory>(parameters_.interaction_parameters, parameters_.use_interactions);
    rate_factory_ = std::make_unique<RateFactory>(input_parameters_.rates);
    cell_factory_ = std::make_unique<CellFactory>(parameters_.site_parameters);
    cell_state_factory_ = std::make_unique<CellStateFactory>(parameters_.site_parameters, rate_factory_.get());
    interaction_state_factory_ = std::make_unique<InteractionStateFactory>(parameters_.interaction_parameters, rate_factory_.get());
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
}

bool CuboidSite::containsPosition(Coordinate3D position) {
    bool contains = (position.x >= lower_bound_.x) && (position.x <= upper_bound_.x) &&
                   (position.y >= lower_bound_.y) && (position.y <= upper_bound_.y) &&
                   (position.z >= lower_bound_.z) && (position.z <= upper_bound_.z);

    return contains;
}

Coordinate3D CuboidSite::getRandomPosition() {

    double   x = random_generator_->generateDouble(lower_bound_.x, upper_bound_.x);
    double   y = random_generator_->generateDouble(lower_bound_.y, upper_bound_.y);
    double   z = random_generator_->generateDouble(lower_bound_.z, upper_bound_.z);
    return Coordinate3D{x, y, z};
}

void CuboidSite::handleBoundaryCross(Agent* agent, Coordinate3D* move_vec, double current_time) {

    boundary_condition_->handleBoundaryCross(agent, move_vec, current_time);
}

Coordinate3D CuboidSite::getLowerLimits() {
    return lower_bound_;
}

Coordinate3D CuboidSite::getUpperLimits() {
    return upper_bound_;
}

Coordinate3D CuboidSite::getRandomBoundaryPoint() {
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

Coordinate3D CuboidSite::generateRandomDirectionVector(Coordinate3D position, double length) {

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

Coordinate3D CuboidSite::generatePersistentDirectionVector(Coordinate3D position,
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

Coordinate3D CuboidSite::generateBackShiftOnContacting(SphereRepresentation* activeSphere, SphereRepresentation* passiveSphere, double mustOverhead) {

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

void CuboidSite::initializeAgents(const abm::util::SimulationParameters::AgentManagerParameters &parameters,
                            double current_time,
                            double time_delta) {

    for (const auto &agent_parameters: parameters.agents) {
        const auto init_distribution = agent_parameters->initial_distribution;
        const auto name = agent_parameters->type;
        input_rates_[name] = 0.0;
        const auto number_of_agents = agent_parameters->number;
        double z_level = (upper_bound_.z + lower_bound_.z)/2;
        for (auto i = 0; i < number_of_agents; ++i) {
            Coordinate3D initial_position = getRandomPosition();
            Coordinate3D initial_vector = Coordinate3D();
            if (dimensions == 2) { initial_position.z = z_level;}

            auto agent = agent_manager_->emplace_back(cell_factory_->createCell(name,
                                                                              std::make_unique<Coordinate3D>(initial_position),
                                                                              agent_manager_->generateNewID(),
                                                                              this,
                                                                              time_delta,
                                                                              current_time, ""));
            if (initial_vector.getMagnitude() != 0) {
                agent->getMovement()->setPreviousMove(&initial_vector);
            }
            if (abm::util::isSubstring("FungalCell", agent->getTypeName())) {
                agent_manager_->addFungalCellToList(agent.get());
            }
        }
    }
}

