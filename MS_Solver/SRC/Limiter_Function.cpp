#include "../INC/Limiter_Function.h"

void MLP_u1_Limiter::set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr)
{
    this->set_of_P1_projected_value_at_vertices_ptr_ = set_of_P1_projected_value_at_vertices_ptr;
}

double MLP_u1_Limiter::limiter_function(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const
{
    const auto criterion_equation_index = criterion.get_criterion_equation_index();
    const auto P1_projected_value_at_vertices = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(cell_index, criterion_equation_index);
    this->set_of_P1_projected_value_at_vertices_ptr_ = &P1_projected_value_at_vertices;
    
    const auto& P0_value_at_vertices = criterion.get_P0_value_at_vertices(cell_index);
    const auto& allowable_min_maxs = criterion.get_criterion_values(cell_index);

    return this->calculate_limiting_value(P0_value_at_vertices, allowable_min_maxs);
};

double MLP_u1_Limiter::limiter_function_opt(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const
{
    const auto& P0_value_at_vertices = criterion.get_P0_value_at_vertices(cell_index);
    const auto& allowable_min_maxs = criterion.get_criterion_values(cell_index);

    return this->calculate_limiting_value(P0_value_at_vertices, allowable_min_maxs);
}

double MLP_u1_Limiter::calculate_limiting_value(const std::vector<double>& P0_value_at_vertices, const std::vector<std::pair<double, double>>& criterion_values) const
{
    const auto num_vertices = this->set_of_P1_projected_value_at_vertices_ptr_->size();

    double limiting_value = 1.0;
    for (ushort j = 0; j < num_vertices; ++j)
    {
        const auto P0_value = P0_value_at_vertices[j];
        const auto [allowable_min, allowable_max] = criterion_values[j];

        const auto P1_mode_criterion_value = this->set_of_P1_projected_value_at_vertices_ptr_->at(j) - P0_value;

        limiting_value = (std::min)(limiting_value, this->calculate_vertex_limiting_value(P0_value, P1_mode_criterion_value, allowable_min, allowable_max));
    }

    return limiting_value;
}

double MLP_u1_Limiter::calculate_vertex_limiting_value(const double P0_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const
{
    if (P1_mode_solution == 0)
    {
        return 1.0;
    }

    if (P1_mode_solution < 0)
    {
        return (std::min)((allowable_min - P0_solution) / P1_mode_solution, 1.0);
    }
    else
    {
        return (std::min)((allowable_max - P0_solution) / P1_mode_solution, 1.0);
    }
};