#pragma once
#include "Indicator.h"

class MLP_u1
{
public:
    MLP_u1(const Grid& gird);

public://Command
    void check(const Discrete_Solution_DG& discrete_solution);    
    
public://Query
    double calculate_limiting_value(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const;

private:
    double calculate_limiting_value(const uint cell_index, const std::vector<double>& P0_value_at_vertices, const std::vector<std::pair<double, double>>& criterion_values) const;
    double calculate_vertex_limiting_value(const double P0_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const;

private:
    ushort criterion_equation_index_ = 0;
    uint num_cells_ = 0;
    std::vector<std::vector<double>> cell_index_to_P1_projected_value_at_vertices_table_;
};
