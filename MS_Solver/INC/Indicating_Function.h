#pragma once
#include "Discrete_Solution.h"

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