#include "../INC/Time_Discrete_Scheme.h"

void SSPRK33::update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const
{  
    auto solution_vw = semi_discrete_equation.solution_vector_wrapper(); //solution vector가 memory reallocation이 될 경우, vetor wrapper의 data ptr들이 dangling pointer가 될 수 있음.
    const auto initial_solution_v = semi_discrete_equation.solution_vector();
    const auto initial_RHS = semi_discrete_equation.calculate_RHS();

    //stage 1
    solution_vw = initial_solution_v + time_step * initial_RHS;
    
    const auto stage1_RHS = semi_discrete_equation.calculate_RHS();

    //stage 2    
    solution_vw = 0.25 * (3 * initial_solution_v + solution_vw + time_step * stage1_RHS);
    
    const auto stage2_RHS = semi_discrete_equation.calculate_RHS();

    //stage 3
    solution_vw = this->c1_3 * (initial_solution_v + 2 * solution_vw + 2 * time_step * stage2_RHS);
}

void SSPRK54::update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const
{
    const auto current_solution_vcw = semi_discrete_equation.solution_vector_constant_wrapper();
    const auto initial_solution_v = semi_discrete_equation.solution_vector();
    const auto initial_RHS = semi_discrete_equation.calculate_RHS();

    //stage 1
    auto stage1_solution_v = initial_solution_v + 0.391752226571890 * time_step * initial_RHS;

    semi_discrete_equation.update_solution(std::move(stage1_solution_v));
    const auto stage1_RHS = semi_discrete_equation.calculate_RHS();

    //stage 2    
    auto stage2_solution_v = 0.444370493651235 * initial_solution_v + 0.555629506348765 * current_solution_vcw + 0.368410593050371 * time_step * stage1_RHS;

    semi_discrete_equation.update_solution(stage2_solution_v);
    const auto stage2_RHS = semi_discrete_equation.calculate_RHS();

    //stage 3
    auto stage3_solution_v = 0.620101851488403 * initial_solution_v + 0.379898148511597 * current_solution_vcw + 0.251891774271694 * time_step * stage2_RHS;

    semi_discrete_equation.update_solution(stage3_solution_v);
    const auto stage3_RHS = semi_discrete_equation.calculate_RHS();

    //stage 4
    auto stage4_solution_v = 0.178079954393132 * initial_solution_v + 0.821920045606868 * current_solution_vcw + 0.544974750228521 * time_step * stage3_RHS;

    semi_discrete_equation.update_solution(std::move(stage4_solution_v));
    const auto stage4_RHS = semi_discrete_equation.calculate_RHS();

    //stage 5
    auto stage5_solution_v = 0.517231671970585 * stage2_solution_v + 0.096059710526147 * stage3_solution_v + 0.063692468666290 * time_step * stage3_RHS + 0.386708617503269 * current_solution_vcw + 0.226007483236906 * time_step * stage4_RHS;

    semi_discrete_equation.update_solution(std::move(stage5_solution_v));
}

std::unique_ptr<Time_Discrete_Scheme> Time_Discrete_Scheme_Factory::make_unique(const Configuration& configuration)
{
    const auto time_discrete_scheme_name = configuration.get("time_discrete_scheme");

    if (ms::contains_icase(time_discrete_scheme_name, "SSPRK", "33"))
    {
        return std::make_unique<SSPRK33>();
    }
    else if (ms::contains_icase(time_discrete_scheme_name, "SSPRK", "54"))
    {
        return std::make_unique<SSPRK54>();
    }
    else
    {
        EXCEPTION("time discrete scheme in configuration is not supproted");
        return nullptr;
    }
};