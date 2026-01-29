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

#ifndef CORE_UTILS_TIMEUTIL_H
#define CORE_UTILS_TIMEUTIL_H

#include <string>
#include <cmath>

namespace abm::util {
    constexpr int standard_fill_width = 13;
    std::string getCurrentLocalTimeAsString();
    std::string generateUniformString(double time, int fill_to_size = standard_fill_width);
}


class SimulationTime {
public:
  // Class for time step management.
    SimulationTime(double delta_t, double max_time) noexcept: max_steps_(
            static_cast<int>((std::ceil(max_time / delta_t)))),
                                                              delta_t_(delta_t),
                                                              last_delta_t_(delta_t),
                                                              max_time_(max_time) {};
    SimulationTime &operator++() noexcept;
    SimulationTime &operator--() noexcept;
    friend int operator%(const SimulationTime &, int);

    void updateTimestep(int new_step) noexcept;
    void updateDeltaT(double new_delta) noexcept;
    [[nodiscard]] double getMaxTime() const noexcept;
    [[nodiscard]] double getCurrentDeltaT() const noexcept;
    [[nodiscard]] double getLastDeltaT() const noexcept;
    [[nodiscard]] double getCurrentTime() const noexcept;
    [[nodiscard]] int getCurrentTimeStep() const noexcept;
    [[nodiscard]] bool endReached() const noexcept;
    [[nodiscard]] bool lastStepBeforEnd() const noexcept;
    [[nodiscard]] int getMaxSteps() const noexcept;
    [[nodiscard]] bool checkForNumberOfExecutions(int number_of_executions, bool execute_first = true) const;
private:
    int time_step_{};
    int max_steps_{};
    double delta_t_{};
    double last_delta_t_{};
    double current_time_{};
    double max_time_{};
};

int operator%(const SimulationTime &, int);

#endif //CORE_UTILS_TIMEUTIL_H
