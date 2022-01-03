#include <cmath>
#include <vector>

#include "Discrete_Solution.h"

class MLP_Criterion_Base
{
public:
    MLP_Criterion_Base(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Query
    ushort get_criterion_equation_index(void) const;
    const std::vector<double>& get_P0_value_at_vertices(const uint cell_index) const;
    const std::vector<std::pair<double, double>>& get_criterion_values(const uint cell_index) const;

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
    std::vector<std::vector<double>> set_of_P0_value_at_vertices;
    std::vector<std::vector<std::pair<double, double>>> set_of_allowable_min_max_criterion_values_;
};

class MLP_Criterion : public MLP_Criterion_Base
{
public:
    MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Command
    void precaclulate(const Discrete_Solution_DG& discrete_solution);

private:
    //precalculate
    std::vector<double> P0_values_;    
};

class Simplex_Decomposed_MLP_Criterion : public MLP_Criterion_Base
{
public:
    Simplex_Decomposed_MLP_Criterion(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Command
    void precaclulate(const Discrete_Solution_DG& discrete_solution);

private:
    std::unordered_map<uint, std::set<uint>> vertex_index_to_matched_vertex_index_set_;

    //precalculate
    std::vector<std::map<uint, double>> set_of_vertex_index_to_simplex_P0_value_;
};

enum class cell_type
{
    normal = 0,
    smooth_extrema = 1,
    trouble = 2,
};

class MLP_Indicator
{
public:
    MLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

public://Command
    void set_precalculated_result(const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr);

public://Query
    cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion& criterion) const;
    cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const Simplex_Decomposed_MLP_Criterion& criterion) const;
    cell_type indicate_opt(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& criterion) const;

private:
    cell_type check_cell_type(const uint cell_index, const double* P1_projected_value_at_vertices, const double* value_at_vertices, const MLP_Criterion_Base& criterion) const;
    bool is_constant(const double value, const double P0_value, const uint cell_index) const;
    bool is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const;
    bool is_smooth_extrema(const double value, const double higher_mode_value, const double P1_mode_value, const double allowable_min, const double allowable_max) const;

private:
    std::vector<double> volumes_;
    std::vector<ushort> set_of_num_vertices_;
        
    //for opt
    const std::vector<double>* set_of_P1_projected_value_at_vertices_ptr_ = nullptr;

    //construction optimization
    static constexpr ushort num_max_vertices = 8;
    mutable std::array<double, num_max_vertices> value_at_vertices_;
    mutable std::array<double, num_max_vertices> P1_projected_value_at_vertices_;
};

class Subcell_Oscillation_Indicator
{
public:
    Subcell_Oscillation_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution);

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

class Shock_Indicator //need to know governing equation!
{
public:
    Shock_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const std::string& governing_equation_name);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution);

public://Query
    bool is_shock(const uint cell_index) const;

private:
    void set_criterion_solution_index(const ushort space_dimension, const std::string& governing_equation_name);

private:    
    short criterion_solution_index_;    
    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;

    //construction optimization
    std::vector<double> average_pressures_;
    std::vector<bool> cell_index_to_is_shock_;    
};

class Discontinuity_Indicator
{
public:
    Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution);

    //test
    const std::vector<double>& get_discontinuity_factor(void) const { return this->discontinuity_factor_; };
    //

private:
    ushort criterion_solution_index_;
    uint num_cells_;

    std::vector<double> cell_index_to_volume_table_;
    std::vector<std::vector<uint>> cell_index_to_face_share_cell_indexes_table_;
    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
    std::vector<bool> cell_index_to_has_discontinuity_table_;

    //for test
    std::vector<double> discontinuity_factor_;
    //

};