#pragma once
#include <vector>

class TIM {};


class SSPRK33 : public TIM {
public:
    template <typename Semi_Discrete_Equation, typename Solution>
    static void update_solutions(std::vector<Solution>& solutions, const double time_step) {
        static const auto num_sol = solutions.size();
        const auto initial_solutions = solutions;

        //stage1
        const auto initial_RHS = Semi_Discrete_Equation::calculate_RHS(initial_solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] += time_step * initial_RHS[i];

        //stage 2
        const auto stage1_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = 0.25 * (3 * initial_solutions[i] + solutions[i] + time_step * stage1_RHS[i]);

        //stage3
        const auto stage2_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = c3_ * (initial_solutions[i] + 2 * solutions[i] + 2 * time_step * stage2_RHS[i]);
    }

private:
    static constexpr double c3_ = 1.0 / 3.0;
};


namespace ms {
    template<typename T>
    inline constexpr bool is_time_integral_method = std::is_base_of_v<TIM, T>;
}