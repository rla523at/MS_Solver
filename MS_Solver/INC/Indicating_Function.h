#pragma once
#include "Stability_Criterion.h"
#include "Measuring_Function.h"

class MLP_Indicating_Function
{
public:
    MLP_Indicating_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution);

public://Query
    Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const;

private:
    Cell_Type check_cell_type(const uint cell_index, const double* value_at_vertices, const MLP_Criterion_Base& criterion) const;
    bool is_constant(const double value, const double P0_value, const uint cell_index) const;
    bool is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const;
    bool is_smooth_extrema(const double value, const double higher_mode_value, const double P1_mode_value, const double allowable_min, const double allowable_max) const;

private:
    ushort criterion_equation_index_;
    uint num_cells_;

    std::vector<double> cell_index_to_volume_reciprocal_table_;
    std::vector<ushort> cell_index_to_num_vertices_table_;
    std::vector<std::vector<double>> cell_index_to_P1_projected_value_at_vertices_table_;

    //construction optimization
    static constexpr ushort num_max_vertices = 8;
    mutable std::array<double, num_max_vertices> value_at_vertices_;
    mutable std::array<double, num_max_vertices> P1_projected_value_at_vertices_;
};

class Subcell_Oscillation_Indicating_Function
{
public:
    Subcell_Oscillation_Indicating_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution);

public://Query
    bool is_typeI_cell(const uint cell_index) const;
    bool is_typeII_cell(const uint cell_index) const;

private:
    ushort criterion_equation_index_ = 0;
    std::vector<ushort> set_of_num_troubled_boundaries_;

    std::vector<double> inner_face_volumes_;
    std::vector<std::pair<uint, uint>> inner_face_oc_nc_index_pairs_;
    std::vector<double> inner_face_characteristic_lengths_;
    std::vector<Euclidean_Vector> set_of_jump_QWs_v_;

    //construction optimization
    static constexpr ushort num_max_jump_QPs = 30;
    mutable std::array<double, num_max_jump_QPs> value_at_ocs_jump_QPs_ = { 0 };
    mutable std::array<double, num_max_jump_QPs> value_at_ncs_jump_QPs_ = { 0 };
    mutable std::array<double, num_max_jump_QPs> value_diff_at_jump_QPs_ = { 0 };
};

class Discontinuity_Indicating_Function
{
public:
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public:
    virtual bool has_discontinuity(const uint cell_index) const abstract;
};