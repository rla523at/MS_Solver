#pragma once
#include "Discrete_Solution.h"

class MLP_Criterion_Base
{
public:
    MLP_Criterion_Base(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);

public://Command
    virtual void precaclulate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    ushort get_criterion_equation_index(void) const;
    const std::vector<double>& get_P0_value_at_vertices(const uint cell_index) const;
    const std::vector<std::pair<double, double>>& get_allowable_min_max_value_at_vertices(const uint cell_index) const;

protected:
    void make_set_of_allowable_min_max_criterion_values(void);

protected:
    ushort criterion_equation_index_ = 0;
    uint num_cells_ = 0;
    const std::unordered_map<uint, std::set<uint>>& vnode_index_to_share_cell_index_set_;
    std::vector<std::vector<uint>> set_of_vertex_indexes_;

    //construction optimization
    static constexpr ushort num_max_vertex_share_cell = 15;
    std::array<double, num_max_vertex_share_cell> criterion_values_;
    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_criterion_value_;
    std::vector<std::vector<double>> cell_index_to_P0_value_at_vertices_table_;
    std::vector<std::vector<std::pair<double, double>>> cell_index_to_allowable_min_max_value_at_vertices_table_;
};