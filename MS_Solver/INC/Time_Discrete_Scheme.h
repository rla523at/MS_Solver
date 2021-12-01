#pragma once
#include "Configuration.h"
#include "Semi_Discrete_Equation.h"

class Time_Discrete_Scheme
{
public:
	virtual void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const abstract;    
};

class SSPRK33 : public Time_Discrete_Scheme
{
public:
    void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const override;

private:
    static constexpr double c1_3 = 1.0 / 3.0;
};

class SSPRK54 : public Time_Discrete_Scheme
{
public:
    void update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const override;

private:
    static constexpr double c1_3 = 1.0 / 3.0;
};

//static class
class Time_Discrete_Scheme_Factory
{
public:
    static std::unique_ptr<Time_Discrete_Scheme> make_unique(const Configuration& configuration);

private:
    Time_Discrete_Scheme_Factory(void) = delete;
};


//class SSPRK54 : public TIM 
//{
//private:
//    SSPRK54(void) = delete;
//
//public:
//    template <typename Semi_Discrete_Equation, typename Solution>
//    static void update_solutions(Semi_Discrete_Equation& semi_discrete_equation, std::vector<Solution>& solutions, const double time_step);
//
//    static std::string name(void) { return "SSPRK54"; };
//};
//
//
//namespace ms {
//    template<typename T>
//    inline constexpr bool is_time_integral_method = std::is_base_of_v<TIM, T>;
//}
//
//
//template <typename Semi_Discrete_Equation, typename Solution>
//static void SSPRK33::update_solutions(Semi_Discrete_Equation& semi_discrete_equation, std::vector<Solution>& solutions, const double time_step) {
    //const auto num_sol = solutions.size();

    //const auto initial_solutions = solutions;
    //const auto initial_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage 1
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] += time_step * initial_RHS[i];
    //semi_discrete_equation.reconstruct(solutions);

    //const auto stage1_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage 2
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = 0.25 * (3 * initial_solutions[i] + solutions[i] + time_step * stage1_RHS[i]);
    //semi_discrete_equation.reconstruct(solutions);
    //
    //const auto stage2_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage 3
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = c1_3 * (initial_solutions[i] + 2 * solutions[i] + 2 * time_step * stage2_RHS[i]);
    //semi_discrete_equation.reconstruct(solutions);
//}
//
//
//template <typename Semi_Discrete_Equation, typename Solution>
//static void SSPRK54::update_solutions(Semi_Discrete_Equation& semi_discrete_equation, std::vector<Solution>& solutions, const double time_step) {
    //const auto num_sol = solutions.size();  

    //const auto initial_solutions = solutions;
    //const auto initial_RHS = semi_discrete_equation.calculate_RHS(solutions);    

    ////stage1
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] += 0.391752226571890 * time_step * initial_RHS[i];    
    //semi_discrete_equation.reconstruct(solutions);    

    //const auto stage1_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage2
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = 0.444370493651235 * initial_solutions[i] + 0.555629506348765 * solutions[i] + 0.368410593050371 * time_step * stage1_RHS[i];
    //semi_discrete_equation.reconstruct(solutions);

    //const auto stage2_solutions = solutions;
    //const auto stage2_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage3
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = 0.620101851488403 * initial_solutions[i] + 0.379898148511597 * solutions[i] + 0.251891774271694 * time_step * stage2_RHS[i];
    //semi_discrete_equation.reconstruct(solutions);

    //const auto stage3_solutions = solutions;
    //const auto stage3_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage4
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = 0.178079954393132 * initial_solutions[i] + 0.821920045606868 * solutions[i] + 0.544974750228521 * time_step * stage3_RHS[i];    
    //semi_discrete_equation.reconstruct(solutions);

    //const auto stage4_RHS = semi_discrete_equation.calculate_RHS(solutions);

    ////stage5
    //for (size_t i = 0; i < num_sol; ++i)
    //    solutions[i] = 0.517231671970585 * stage2_solutions[i] + 0.096059710526147 * stage3_solutions[i] + 0.063692468666290 * time_step * stage3_RHS[i] + 0.386708617503269 * solutions[i] + 0.226007483236906 * time_step * stage4_RHS[i];
    //semi_discrete_equation.reconstruct(solutions);
//}