#pragma once
#include "Measuring_Function.h"

class Scaled_Difference : public Discontinuity_Measuring_Function
{
public:
    Scaled_Difference(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    std::vector<double> measure(const Discrete_Solution_DG& discrete_solution) const override;

private:
    ushort criterion_solution_index_;
    uint num_cells_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
};

class Extrapolation_Difference : public Discontinuity_Measuring_Function
{
public:
    Extrapolation_Difference(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    std::vector<double> measure(const Discrete_Solution_DG& discrete_solution) const override;    

private:
    ushort criterion_solution_index_;
    uint num_cells_;
    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
    std::vector<double> cell_index_to_volume_reciprocal_table_;
};

class Difference_Of_Extrapolatiion_Difference : public Discontinuity_Measuring_Function
{
public:
    Difference_Of_Extrapolatiion_Difference(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
        : num_cells_(grid.num_cells())
        , extrapolate_difference_(grid, discrete_solution, criterion_solution_index)
        , cell_index_to_face_share_cell_indexes_table_(grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry()) {};

public:
    std::vector<double> measure(const Discrete_Solution_DG& discrete_solution) const override;

private:
    uint num_cells_ = 0;
    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
    Extrapolation_Difference extrapolate_difference_;
};

class Max_Divergence_of_Velocity : public Discontinuity_Measuring_Function
{
public:
    Max_Divergence_of_Velocity(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    std::vector<double> measure(const Discrete_Solution_DG& discrete_solution) const override;

private:
    uint num_cells_ = 0;
    std::vector<ushort> cell_index_to_num_QPs_;

    //construction
    static constexpr int max_QPs = 30;
    static inline std::array<Euclidean_Vector, max_QPs> solution_at_cell_QPs;
    static inline std::array<Euclidean_Vector, max_QPs> ddx_GE_solution_at_cell_QPs;
    static inline std::array<Euclidean_Vector, max_QPs> ddy_GE_solution_at_cell_QPs;
};