#pragma once
#include "Indicator.h"

class MLP_u1_Limiter
{
public://Command
    void set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr);
    
public://Query
    double limiter_function(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const;
    double limiter_function_opt(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const;

private:
    double calculate_limiting_value(const std::vector<double>& P0_value_at_vertices, const std::vector<std::pair<double, double>>& criterion_values) const;
    double calculate_vertex_limiting_value(const double P0_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const;

private:
    //for optimization
    mutable const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr_ = nullptr;
};
