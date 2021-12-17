#include <cmath>
#include <vector>

#include "Discrete_Solution.h"

enum class cell_type
{
    normal,
    smooth_extrema,
    trouble
};

class MLP_Criterion
{
public:
    MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Command
    void caclulate(const Discrete_Solution_DG& discrete_solution);

public://Query
    const std::vector<std::pair<double, double>>& get_criterion_values(const uint cell_index) const;
    double get_P0_value(const uint cell_index) const;
    ushort get_criterion_equation_index(void) const;

private:
    ushort criterion_equation_index_ = 0;
    uint num_cells_ = 0;
    const std::unordered_map<uint, std::set<uint>>& vnode_index_to_share_cell_index_set_;
    std::vector<std::vector<uint>> set_of_vnode_indexes_;

    //precalculate
    std::vector<double> P0_values_;
    std::vector<std::vector<std::pair<double, double>>> set_of_allowable_min_max_criterion_values_;

    //construction optimization
    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value_;
};

class MLP_Indicator
{
public:
    MLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Command
    void set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr);

public://Query
    cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const;
    cell_type indicate_opt(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const;

private:
    cell_type check_cell_type(const uint cell_index, const double* P1_projected_value_at_vertices, const double* value_at_vertices, const MLP_Criterion& criterion) const;
    bool is_constant(const double solution, const double P0_solution, const uint cell_index) const;
    bool is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const;
    bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) const;

private:
    std::vector<double> volumes_;
    std::vector<ushort> set_of_num_vertices_;
        
    //for optimization
    const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr_ = nullptr;
};