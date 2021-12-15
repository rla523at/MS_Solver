#include "../INC/Time_Discrete_Scheme.h"

void SSPRK33::update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const
{  
    constexpr auto c1_3 = 1.0 / 3.0;

    auto solution_vw = semi_discrete_equation.discrete_solution_vector_wrapper();
    const auto initial_solution_v = semi_discrete_equation.discrete_solution_vector();
    const auto initial_RHS = semi_discrete_equation.calculate_RHS();

    //stage 1
    solution_vw = initial_solution_v + time_step * initial_RHS;
    
    const auto stage1_RHS = semi_discrete_equation.calculate_RHS();

    //stage 2    
    solution_vw = 0.25 * (3 * initial_solution_v + solution_vw + time_step * stage1_RHS);
    
    const auto stage2_RHS = semi_discrete_equation.calculate_RHS();

    //stage 3
    solution_vw = c1_3 * (initial_solution_v + 2 * solution_vw + 2 * time_step * stage2_RHS);
}

void SSPRK54::update(Semi_Discrete_Equation& semi_discrete_equation, const double time_step) const
{
    auto solution_vw = semi_discrete_equation.discrete_solution_vector_wrapper();
    const auto initial_solution_v = semi_discrete_equation.discrete_solution_vector();
    const auto initial_RHS = semi_discrete_equation.calculate_RHS();

    //stage 1
    solution_vw = initial_solution_v + 0.391752226571890 * time_step * initial_RHS;

    semi_discrete_equation.reconstruct();
    const auto stage1_RHS = semi_discrete_equation.calculate_RHS();

    //stage 2    
    solution_vw = 0.444370493651235 * initial_solution_v + 0.555629506348765 * solution_vw + 0.368410593050371 * time_step * stage1_RHS;

    semi_discrete_equation.reconstruct();
    const auto stage2_solution_v = semi_discrete_equation.discrete_solution_vector();
    const auto stage2_RHS = semi_discrete_equation.calculate_RHS();

    //stage 3
    solution_vw = 0.620101851488403 * initial_solution_v + 0.379898148511597 * solution_vw + 0.251891774271694 * time_step * stage2_RHS;

    semi_discrete_equation.reconstruct();
    const auto stage3_solution_v = semi_discrete_equation.discrete_solution_vector();
    const auto stage3_RHS = semi_discrete_equation.calculate_RHS();

    //stage 4
    solution_vw = 0.178079954393132 * initial_solution_v + 0.821920045606868 * solution_vw + 0.544974750228521 * time_step * stage3_RHS;

    semi_discrete_equation.reconstruct();
    const auto stage4_RHS = semi_discrete_equation.calculate_RHS();

    //stage 5
    solution_vw = 0.517231671970585 * stage2_solution_v + 0.096059710526147 * stage3_solution_v + 0.063692468666290 * time_step * stage3_RHS + 0.386708617503269 * solution_vw + 0.226007483236906 * time_step * stage4_RHS;

    semi_discrete_equation.reconstruct();
}

std::unique_ptr<Time_Discrete_Scheme> Time_Discrete_Scheme_Factory::make_unique(const Configuration& configuration)
{
    const auto& time_discrete_scheme_name = configuration.get_time_discrete_scheme();

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