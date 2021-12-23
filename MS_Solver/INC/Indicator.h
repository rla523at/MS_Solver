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
    std::vector<std::vector<uint>> set_of_face_share_cell_indexes_;

    //construction optimization
    std::vector<double> average_pressures_;
    std::vector<bool> are_shock_;    
};

class Discontinuity_Indicator
{
public:
    Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
        : criterion_solution_index_(criterion_solution_index)
    {
        this->num_cells_ = grid.num_cells();
        
        this->has_discontinuities_.resize(this->num_cells_, false);
        this->set_of_QWs_v_.resize(this->num_cells_);
        this->set_of_face_share_cell_indexes_ = grid.set_of_face_share_cell_indexes_consider_pbdry();

        
        std::vector<Quadrature_Rule> quadrature_rules(this->num_cells_);

        for (uint i = 0; i < this->num_cells_; ++i)
        {
            const auto solution_degree = discrete_solution.solution_degree(i);
            quadrature_rules[i] = grid.get_cell_quadrature_rule(i, solution_degree);
            this->set_of_QWs_v_[i] = quadrature_rules[i].weights;
        }

        //for precalculation
        discrete_solution.precalculate_set_of_cell_index_to_target_cell_basis_QPs_m_(quadrature_rules, this->set_of_face_share_cell_indexes_);
        //
    }

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution)
    {
        static constexpr ushort max_num_QPs = 100;
        std::array<double, max_num_QPs> nth_solution_at_target_cell_QPs = { 0 };

        //for test
        std::vector<double> discontinuity_factor_(this->num_cells_);


        for (uint i = 0; i < this->num_cells_; ++i)
        {
            const auto P0_value = discrete_solution.calculate_P0_nth_solution(i, this->criterion_solution_index_);

            const auto& QWs_v = this->set_of_QWs_v_[i];
            const auto num_QPs = static_cast<int>(QWs_v.size());

            const auto& face_share_cell_indexes = this->set_of_face_share_cell_indexes_[i];
            const auto num_face_share_cells = face_share_cell_indexes.size();


            for (ushort j = 0; j < num_face_share_cells; ++j)
            {
                const auto my_cell_index = face_share_cell_indexes[j];

                discrete_solution.calculate_nth_solution_at_target_cell_QPs(nth_solution_at_target_cell_QPs.data(), i, my_cell_index, this->criterion_solution_index_);
                
                const auto P0_value_by_extrapolate = ms::BLAS::x_dot_y(num_QPs, nth_solution_at_target_cell_QPs.data(), QWs_v.data());

                discontinuity_factor_[i] += std::abs(P0_value_by_extrapolate - P0_value);
            }

            discontinuity_factor_[i] /= num_face_share_cells;
        }
    }

private:
    ushort criterion_solution_index_;
    uint num_cells_;

    std::vector<std::vector<uint>> set_of_face_share_cell_indexes_;
    std::vector<Euclidean_Vector> set_of_QWs_v_;
    std::vector<bool> has_discontinuities_;
};