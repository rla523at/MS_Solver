#include "../INC/Time_Discrete_Scheme.h"

void SSPRK33::update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const
{    
    const auto& current_solution_vcw = semi_discrete_equation.solution_vector_constant_wrapper();
    const auto initial_solution_v = semi_discrete_equation.solution_vector();
    const auto initial_RHS = semi_discrete_equation.calculate_RHS();

    //stage 1
    auto stage1_solution_v = initial_solution_v + time_step * initial_RHS;

    semi_discrete_equation.update_solution(std::move(stage1_solution_v));
    const auto stage1_RHS = semi_discrete_equation.calculate_RHS();

    //stage 2    
    auto stage2_solution_v = 0.25 * (3 * initial_solution_v + current_solution_vcw + time_step * stage1_RHS);
    
    semi_discrete_equation.update_solution(std::move(stage2_solution_v));
    const auto stage2_RHS = semi_discrete_equation.calculate_RHS();

    //stage 3
    auto stage3_solution_v = this->c1_3 * (initial_solution_v + 2 * current_solution_vcw + 2 * time_step * stage2_RHS);
    semi_discrete_equation.update_solution(std::move(stage3_solution_v));
}

std::unique_ptr<Time_Discrete_Scheme> Time_Discrete_Scheme_Factory::make_unique(const Configuration& configuration)
{
    const auto time_discrete_scheme_name = configuration.get("time_discrete_scheme");

    if (ms::contains_icase(time_discrete_scheme_name, "SSPRK", "33"))
    {
        return std::make_unique<SSPRK33>();
    }
    else
    {
        EXCEPTION("time discrete scheme in configuration is not supproted");
        return nullptr;
    }
};