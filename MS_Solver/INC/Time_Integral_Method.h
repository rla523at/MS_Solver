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


class SSPRK54 : public TIM {
public:
    template <typename Semi_Discrete_Equation, typename Solution>
    static void update_solutions(std::vector<Solution>& solutions, const double time_step) {
        static const auto num_sol = solutions.size();

        //auto stage1_solutions = solutions;
        //auto stage2_solutions = solutions;
        //auto stage3_solutions = solutions;
        //auto stage4_solutions = solutions;


        //const auto initial_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        ////stage1
        //for (size_t i = 0; i < num_sol; ++i)
        //    stage1_solutions[i] = solutions[i] + 0.391752226571890 * time_step * initial_RHS[i];
        //const auto stage1_RHS = Semi_Discrete_Equation::calculate_RHS(stage1_solutions);

        ////stage2
        //for (size_t i = 0; i < num_sol; ++i)
        //    stage2_solutions[i] = 0.444370493651235 * solutions[i] + 0.555629506348765 * stage1_solutions[i] + 0.368410593050371 * time_step * stage1_RHS[i];
        //const auto stage2_RHS = Semi_Discrete_Equation::calculate_RHS(stage2_solutions);

        ////stage3
        //for (size_t i = 0; i < num_sol; ++i)
        //    stage3_solutions[i] = 0.620101851488403 * solutions[i] + 0.379898148511597 * stage2_solutions[i] + 0.251891774271694 * time_step * stage2_RHS[i];
        //const auto stage3_RHS = Semi_Discrete_Equation::calculate_RHS(stage3_solutions);

        ////stage4
        //for (size_t i = 0; i < num_sol; ++i)
        //    stage4_solutions[i] = 0.178079954393132 * solutions[i] + 0.821920045606868 * stage3_solutions[i] + 0.544974750228521 * time_step * stage3_RHS[i];
        //const auto stage4_RHS = Semi_Discrete_Equation::calculate_RHS(stage4_solutions);

        ////stage5
        //for (size_t i = 0; i < num_sol; ++i)
        //    solutions[i] = 0.517231671970585 * stage2_solutions[i] + 0.096059710526147 * stage3_solutions[i] + 0.063692468666290 * time_step * stage3_RHS[i] + 0.386708617503269 * stage4_solutions[i] + 0.226007483236906 * time_step * stage4_RHS[i];


        const auto initial_solutions = solutions;

        //stage1
        const auto initial_RHS = Semi_Discrete_Equation::calculate_RHS(initial_solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] += 0.391752226571890 * time_step * initial_RHS[i];

        //stage2
        const auto stage1_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = 0.444370493651235 * initial_solutions[i] + 0.555629506348765 * solutions[i] + 0.368410593050371 * time_step * stage1_RHS[i];

        const auto stage2_solutions = solutions;

        //stage3
        const auto stage2_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = 0.620101851488403 * initial_solutions[i] + 0.379898148511597 * solutions[i] + 0.251891774271694 * time_step * stage2_RHS[i];

        const auto stage3_solutions = solutions;

        //stage4
        const auto stage3_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = 0.178079954393132 * initial_solutions[i] + 0.821920045606868 * solutions[i] + 0.544974750228521 * time_step * stage3_RHS[i];

        //stage5
        const auto stage4_RHS = Semi_Discrete_Equation::calculate_RHS(solutions);
        for (size_t i = 0; i < num_sol; ++i)
            solutions[i] = 0.517231671970585 * stage2_solutions[i] + 0.096059710526147 * stage3_solutions[i] + 0.063692468666290 * time_step * stage3_RHS[i] + 0.386708617503269 * solutions[i] + 0.226007483236906 * time_step * stage4_RHS[i];
    }
};


namespace ms {
    template<typename T>
    inline constexpr bool is_time_integral_method = std::is_base_of_v<TIM, T>;
}