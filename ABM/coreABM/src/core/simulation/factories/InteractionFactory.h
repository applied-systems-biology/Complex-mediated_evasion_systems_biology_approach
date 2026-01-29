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

#ifndef CORE_SIMULATION_INTERACTIONFACTORY_H
#define CORE_SIMULATION_INTERACTIONFACTORY_H

#include <memory>
#include <utility>
#include <map>

#include "core/simulation/Interaction.h"
#include "core/analyser/Analyser.h"

class Collision;
class Cell;

class InteractionFactory {
public:
  // Factory class for the interactions specified in the simulerator configuration
    InteractionFactory(
            const std::vector<std::unique_ptr<abm::util::SimulationParameters::InteractionParameters>> &interaction_parameters, bool use_interactions);

    unsigned int generateInteractionId();

    virtual std::shared_ptr<Interaction> createInteraction(double time_delta,
                                                          double current_time,
                                                          const std::shared_ptr<Collision> &collision,
                                                          InSituMeasurements *measurements);
    std::shared_ptr<Interaction> createAvoidanceInteraction(Cell *cell_1,
                                                                   std::shared_ptr<Collision> collision,
                                                                   double time_delta,
                                                                   double current_time);
    bool isInteractionsOn();

private:
    unsigned int interaction_id_;
    std::map<std::string, std::string> interaction_types_;
    std::map<std::string, std::vector<std::pair<std::string, std::vector<std::string>>>> interaction_conditions_;
    std::map<std::pair<std::string, std::string>, std::string> interaction_pair_types_;

  protected:
    std::tuple<std::string, std::string> retrieveInteractionIdentifier(Cell *cell_1, Cell *cell_2);
};

#endif /* CORE_SIMULATION_INTERACTIONFACTORY_H */
