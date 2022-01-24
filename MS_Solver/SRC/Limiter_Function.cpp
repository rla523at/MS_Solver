#include "../INC/Limiting_Function.h"

MLP_u1::MLP_u1(const Grid& grid)
{
    this->num_cells_ = grid.num_cells();
    this->cell_index_to_P1_projected_value_at_vertices_table_.resize(this->num_cells_);

    const auto cell_index_to_num_vertices_table = grid.cell_index_to_num_vertices_table();
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index].resize(cell_index_to_num_vertices_table[cell_index]);
    }
}

void MLP_u1::check(const Discrete_Solution_DG& discrete_solution)
{
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        discrete_solution.calculate_P1_projected_nth_solution_at_vertices(this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index].data(), cell_index, this->criterion_equation_index_);
    }
}

//double MLP_u1_Limiter::limiter_function(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const
//{
//    const auto criterion_equation_index = criterion.get_criterion_equation_index();
//    const auto P1_projected_value_at_vertices = discrete_solution.calculate_P1_projected_nth_solution_at_vertices(cell_index, criterion_equation_index);
//    this->set_of_P1_projected_value_at_vertices_ptr_ = &P1_projected_value_at_vertices;
//    
//    const auto& P0_value_at_vertices = criterion.get_P0_value_at_vertices(cell_index);
//    const auto& allowable_min_maxs = criterion.get_criterion_values(cell_index);
//
//    return this->calculate_limiting_value(P0_value_at_vertices, allowable_min_maxs);
//};

double MLP_u1::calculate_limiting_value(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const
{
    const auto& P0_value_at_vertices = criterion.get_P0_value_at_vertices(cell_index);
    const auto& allowable_min_maxs = criterion.get_allowable_min_max_value_at_vertices(cell_index);

    return this->calculate_limiting_value(cell_index, P0_value_at_vertices, allowable_min_maxs);
}

double MLP_u1::calculate_limiting_value(const uint cell_index, const std::vector<double>& P0_value_at_vertices, const std::vector<std::pair<double, double>>& criterion_values) const
{
    const auto& P1_projected_value_at_vertices = this->cell_index_to_P1_projected_value_at_vertices_table_[cell_index];
    const auto num_vertices = P1_projected_value_at_vertices.size();

    double limiting_value = 1.0;
    for (ushort j = 0; j < num_vertices; ++j)
    {
        const auto P0_value = P0_value_at_vertices[j];
        const auto [allowable_min, allowable_max] = criterion_values[j];

        const auto P1_mode_criterion_value = P1_projected_value_at_vertices[j] - P0_value;

        limiting_value = (std::min)(limiting_value, this->calculate_vertex_limiting_value(P0_value, P1_mode_criterion_value, allowable_min, allowable_max));
    }

    return limiting_value;
}

double MLP_u1::calculate_vertex_limiting_value(const double P0_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const
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