#pragma once
#include "Discrete_Solution.h"

class Scaled_Average_Difference_Measuring_Function
{
public:
    Scaled_Average_Difference_Measuring_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    std::vector<double> measure_infc_index_to_scaled_average_difference_table(const Discrete_Solution_DG& discrete_solution) const;

private:
    ushort criterion_solution_index_;
    uint num_infcs_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
};

class Extrapolation_Differences_Measuring_Function
{
public:
    Extrapolation_Differences_Measuring_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    std::vector<std::vector<double>> measure_cell_index_to_extrapolation_differences(const Discrete_Solution_DG& discrete_solution) const;    

private:
    ushort criterion_solution_index_;
    uint num_cells_;
    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
    std::vector<double> cell_index_to_volume_reciprocal_table_;
};

class Divergence_Velocity_Measuring_Function
{
public:
    Divergence_Velocity_Measuring_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public:
    std::vector<std::vector<double>> measure_cell_index_to_divergence_velocities_table(const Discrete_Solution_DG& discrete_solution) const;

private:
    uint num_cells_ = 0;
    std::vector<ushort> cell_index_to_num_QPs_;

    //construction optimization
    static constexpr int max_QPs = 100;
    static inline std::array<Euclidean_Vector, max_QPs> solution_at_cell_QPs;
    static inline std::array<Euclidean_Vector, max_QPs> ddx_GE_solution_at_cell_QPs;
    static inline std::array<Euclidean_Vector, max_QPs> ddy_GE_solution_at_cell_QPs;
};

class Average_Solution_Jump_Measuring_Function
{
public:
    Average_Solution_Jump_Measuring_Function(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public://Command
    std::vector<double> measure_infc_index_to_average_solution_jump_table(const Discrete_Solution_DG& discrete_solution);
    
private:
    ushort criterion_solution_index_ = 0;
    uint num_infcs_ = 0;

    std::vector<double> infc_index_to_reciprocal_volume_table_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;
    std::vector<Euclidean_Vector> infc_index_to_jump_QWs_v_table_;

    //construction optimization
    static constexpr ushort num_max_jump_QPs = 30;
    mutable std::array<double, num_max_jump_QPs> value_at_ocs_jump_QPs_ = { 0 };
    mutable std::array<double, num_max_jump_QPs> value_at_ncs_jump_QPs_ = { 0 };
    mutable std::array<double, num_max_jump_QPs> value_diff_at_jump_QPs_ = { 0 };
};




//class Difference_Of_Extrapolatiion_Difference
//{
//public:
//    Difference_Of_Extrapolatiion_Difference(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
//        : num_cells_(grid.num_cells())
//        , extrapolate_difference_(grid, discrete_solution, criterion_solution_index)
//        , cell_index_to_face_share_cell_indexes_table_(grid.cell_index_to_face_share_cell_indexes_table_consider_pbdry()) {};
//
//public:
//    std::vector<double> measure_infc_index_to_scaled_average_difference_table(const Discrete_Solution_DG& discrete_solution) const;
//
//private:
//    uint num_cells_ = 0;
//    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
//    Extrapolation_Differences_Measuring_Function extrapolate_difference_;
//};