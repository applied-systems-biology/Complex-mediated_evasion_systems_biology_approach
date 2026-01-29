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


#ifndef CORE_SIMULATION_FungalCellExample_H
#define CORE_SIMULATION_FungalCellExample_H

#include "core/simulation/cells/FungalCell.h"
#include "io_utils_example.h"


class FungalCellExample : public FungalCell {
public:
    // To create your own Celltypes, you have to create your own CellFactory, to make sure the cell can be created
    // This must be included in the constructor of the Site (here: CuboidSiteExample)
    // The corresponding CellType must be used in the simulator-config.json etc.
    FungalCellExample(std::unique_ptr<Coordinate3D> c, int id, Site *site,
               double time_delta, double current_time) : FungalCell(std::move(c), id, site, time_delta, current_time) {}

    void doMorphologicalChanges(double timestep, double current_time) final;

    // To guarantee the same functionality for derived class "FungalCellExample" as for "FungalCell":
    // Make sure the getTypeName returns a substring of "FungalCell"
    std::string getTypeName() { return "FungalCellExample"; };

    void setup(double time_delta, double current_time, abm::util::SimulationParameters::AgentParameters* parameters) final;
    abm::utilExample::FungalParametersExample* getFungalParameter(){return fc_parameters;};

private:
    abm::utilExample::FungalParametersExample* fc_parameters;
};

#endif /* CORE_SIMULATION_FungalCellExample_H */
