#pragma once
#include "Reconstruction_Method_FVM.h"
#include "Polynomial.h"

class HOM_Reconstruction : public RM {};

template <ushort space_dimension_, ushort solution_order_>
class Polynomial_Reconstruction : public HOM_Reconstruction
{
protected:
    using This_ = Polynomial_Reconstruction<space_dimension_, solution_order_>;
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;

    static constexpr ushort num_basis_ = ms::combination_with_repetition(1 + space_dimension_, solution_order_);

protected:
    inline static std::vector<Vector_Function<Polynomial<space_dimension_>, num_basis_>> set_of_basis_functions_;

private:
    Polynomial_Reconstruction(void) = delete;

public:
    static void initialize(const Grid<space_dimension_>& grid);
    static auto calculate_set_of_transposed_gradient_basis(void);
    static auto calculate_basis_node(const uint cell_index, const Space_Vector_& node);
    static Dynamic_Matrix calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes);
    static std::vector<Dynamic_Matrix> calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes);
    static double calculate_P0_basis_value(const uint cell_index, const Space_Vector_& center_node);
    static std::vector<double> calculate_P0_basis_values(const std::vector<Space_Vector_>& center_nodes);
    static constexpr ushort num_basis(void) { return This_::num_basis_; };
    static constexpr ushort solution_order(void) { return solution_order_; };
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static std::string name(void) { return "P" + std::to_string(solution_order_) + "_Polynomial_Reconstruction"; };
};


template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
class hMLP_Reconstruction : public Polynomial_Reconstruction<space_dimension_, solution_order_>
{
private:
    static_require(1 <= solution_order_, "solution order should be greater than 1 for hMLP reconstruction");

    using Base_ = Polynomial_Reconstruction<space_dimension_, solution_order_>;
    using This_ = hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>;
    using Solution_Coefficients_ = Matrix<num_equation, This_::num_basis_>;

    static constexpr ushort criterion_variable_index_ = 0; //rho

private:
    inline static std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_indexes_;
    inline static std::vector<std::vector<uint>> set_of_vnode_indexes_;
    inline static std::vector<Dynamic_Matrix> P1_mode_basis_vnodes_;
    inline static std::vector<Dynamic_Matrix> P1_projected_basis_vnodes_;
    inline static std::vector<double> P0_basis_values_;
    inline static std::vector<Dynamic_Matrix> basis_vnodes_;
    inline static std::array<ushort, solution_order_ + 1> num_Pn_basis_;

private:
    hMLP_Reconstruction(void) = delete;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    static void reconstruct(std::vector<Solution_Coefficients_>& solutions);
    static std::string name(void) { return "P" + std::to_string(solution_order_) + "_hMLP_Reconstruction"; };

    //private: //for test
    static auto calculate_vertex_node_index_to_allowable_min_max_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions);
    static bool is_satisfy_P1_projected_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max);
    static bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max);
};


//template definition part
template <ushort space_dimension_, ushort solution_order_>
void Polynomial_Reconstruction<space_dimension_, solution_order_>::initialize(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    This_::set_of_basis_functions_.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        This_::set_of_basis_functions_.push_back(cell_geometry.orthonormal_basis_functions<solution_order_>());
    }

    Log::content_ << std::left << std::setw(50) << "@ Polynomial reconstruction precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_transposed_gradient_basis(void) {
    const auto num_cell = This_::set_of_basis_functions_.size();

    std::vector<Matrix_Function<Polynomial<space_dimension_>, space_dimension_, This_::num_basis_>> transposed_basis_Jacobians(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        auto& transposed_jacobian_basis = transposed_basis_Jacobians[i];
        const auto& basis_function = This_::set_of_basis_functions_[i];

        for (ushort j = 0; j < This_::num_basis_; ++j)
            transposed_jacobian_basis.change_column(j, basis_function[j].gradient());
    }

    return transposed_basis_Jacobians;
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_node(const uint cell_index, const Space_Vector_& node) {
    const auto& basis_functions = This_::set_of_basis_functions_[cell_index];
    return basis_functions(node);
}


template <ushort space_dimension_, ushort solution_order_>
Dynamic_Matrix Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) {
    const auto& basis_functions = This_::set_of_basis_functions_[cell_index];

    const auto num_node = nodes.size();

    Dynamic_Matrix basis_nodes(This_::num_basis_, num_node);
    for (size_t j = 0; j < num_node; ++j)
        basis_nodes.change_column(j, basis_functions(nodes[j]));

    return basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
double Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_value(const uint cell_index, const Space_Vector_& center_node) {
    const auto& basis_function = This_::set_of_basis_functions_[cell_index];
    const auto& P0_basis_function = basis_function[0];

    return P0_basis_function(center_node);
}

template <ushort space_dimension_, ushort solution_order_>
std::vector<Dynamic_Matrix> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes) {
    const auto num_cell = This_::set_of_basis_functions_.size();

    dynamic_require(set_of_nodes.size() == num_cell, "set size should be same with num cell");

    std::vector<Dynamic_Matrix> set_of_basis_nodes;
    set_of_basis_nodes.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i)
        set_of_basis_nodes.push_back(This_::calculate_basis_nodes(i, set_of_nodes[i]));


    return set_of_basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
std::vector<double> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_values(const std::vector<Space_Vector_>& center_nodes) {
    const auto num_cell = center_nodes.size();
    std::vector<double> P0_basis_values(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& basis_function = This_::set_of_basis_functions_[i];
        const auto& P0_basis_function = basis_function[0];
        const auto& cell_center = center_nodes[i];

        P0_basis_values[i] = P0_basis_function(cell_center);
    }

    return P0_basis_values;
}

template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
void hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::initialize(Grid<space_dimension_>&& grid) {    
    SET_TIME_POINT;
    
    for (ushort i = 0; i <= solution_order_; ++i) 
        This_::num_Pn_basis_[i] = ms::combination_with_repetition(1 + space_dimension_, i);

    Base_::initialize(grid);
    This_::vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_indexes);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto num_cell = cell_elements.size();
    This_::set_of_vnode_indexes_.reserve(num_cell);
    This_::basis_vnodes_.reserve(num_cell);
    This_::P1_projected_basis_vnodes_.reserve(num_cell);
    This_::P1_mode_basis_vnodes_.reserve(num_cell);
    This_::P0_basis_values_.reserve(num_cell);

    constexpr ushort P1_solution_order = 1;
    constexpr auto num_P1_projected_basis = ms::combination_with_repetition(1 + space_dimension_, P1_solution_order);
    constexpr auto num_over_P2_basis = This_::num_basis_ - num_P1_projected_basis;

    constexpr ushort P0_basis_row_index = 0;

    for (uint i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        This_::set_of_vnode_indexes_.push_back(cell_element.vertex_node_indexes());

        const auto vnodes = cell_geometry.vertex_nodes();
        const auto num_vnode = vnodes.size();

        auto basis_vnode = This_::calculate_basis_nodes(i, vnodes);

        auto P1_projected_basis_vnode = basis_vnode;
        if (solution_order_ > 1) {
            Dynamic_Matrix P1_projection_matrix(num_over_P2_basis, num_vnode);
            P1_projected_basis_vnode.change_rows(num_P1_projected_basis, P1_projection_matrix);
        }

        Dynamic_Matrix P1_mode_basis_vnode = P1_projected_basis_vnode;
        P1_mode_basis_vnode.change_row(P0_basis_row_index, Dynamic_Euclidean_Vector(num_vnode));

        const auto P0_basis_value = P1_projected_basis_vnode.at(P0_basis_row_index, 0);

        This_::basis_vnodes_.push_back(std::move(basis_vnode));
        This_::P1_projected_basis_vnodes_.push_back(std::move(P1_projected_basis_vnode));
        This_::P1_mode_basis_vnodes_.push_back(std::move(P1_mode_basis_vnode));
        This_::P0_basis_values_.push_back(P0_basis_value);
    }

    Log::content_ << std::left << std::setw(50) << "@ hMLP precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
void hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::reconstruct(std::vector<Solution_Coefficients_>& solution_coefficients) {
    const auto num_cell = solution_coefficients.size();

    std::vector<Euclidean_Vector<num_equation>> P0_solutions(num_cell);
    for (uint i = 0; i < num_cell; ++i)
        P0_solutions[i] = (solution_coefficients[i] * This_::P0_basis_values_[i]).column(0);

    const auto vnode_index_to_allowable_min_max_solution = This_::calculate_vertex_node_index_to_allowable_min_max_solution(P0_solutions);

    for (uint i = 0; i < num_cell; ++i) {
        auto temporal_solution_order = solution_order_;

        auto& solution_coefficient = solution_coefficients[i];

        auto solution_vnode = solution_coefficient * This_::basis_vnodes_[i];
        auto P1_projected_solution_vnode = solution_coefficient * This_::P1_projected_basis_vnodes_[i];

        const auto& vnode_indexes = This_::set_of_vnode_indexes_[i];
        const auto num_vnode = vnode_indexes.size();

        bool end_limiting = false;

        for (ushort j = 0; j < num_vnode; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_solution.at(vnode_index);

            while (true) {
                const auto P1_projected_solution = P1_projected_solution_vnode.column<num_equation>(j);
                const auto P1_projected_criterion_variable = P1_projected_solution[This_::criterion_variable_index_];

                if (This_::is_satisfy_P1_projected_MLP_condition(P1_projected_criterion_variable, allowable_min, allowable_max))                     
                    break;                    
                else {
                    const auto solution = solution_vnode.column<num_equation>(j);
                    const auto solution_criterion_variable = solution[This_::criterion_variable_index_];

                    const auto higher_mode_criterion_variable = solution_criterion_variable - P1_projected_criterion_variable;
                    const auto P1_mode_criterion_variable = P1_projected_criterion_variable - P0_solutions[i][This_::criterion_variable_index_];
                    const auto P0_mode_criterion_variable = P0_solutions[i][This_::criterion_variable_index_];

                    if (This_::is_smooth_extrema(solution_criterion_variable, higher_mode_criterion_variable, P1_mode_criterion_variable, allowable_min, allowable_max))
                        break;

                    else {
                        if (temporal_solution_order == 1) {
                            const auto limiting_value = MLP_u1_Limiting_Strategy::limit(P1_mode_criterion_variable, P0_mode_criterion_variable, allowable_min, allowable_max);
                            
                            std::array<double, This_::num_basis_> limiting_values;
                            limiting_values.fill(limiting_value);
                            limiting_values[0] = 1.0; // keep P0
                            
                            const auto limiting_matrix = Matrix<This_::num_basis_, This_::num_basis_>::diagonal_matrix(limiting_values);

                            solution_coefficient *= limiting_matrix;
                            end_limiting = true;
                            break;
                        }
                        else {
                            //limiting high order term
                            const auto num_Pnm1_basis = This_::num_Pn_basis_[temporal_solution_order - 1];
                            const auto num_Pn_basis = This_::num_Pn_basis_[temporal_solution_order];
                            const auto num_delete_basis = num_Pn_basis - num_Pnm1_basis;
                            const Dynamic_Matrix order_down_matrix(num_equation, num_delete_basis);
                            solution_coefficient.change_columns(num_Pnm1_basis, order_down_matrix);     

                            //re-calculate limited solution
                            solution_vnode = solution_coefficient * This_::basis_vnodes_[i];
                            P1_projected_solution_vnode = solution_coefficient * This_::P1_projected_basis_vnodes_[i];
                        }
                    }
                }    
            }

            if (end_limiting == true)
                break;
        }
    }
}


template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
auto hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::calculate_vertex_node_index_to_allowable_min_max_solution(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) {
    const auto num_solution = P0_solutions.size();
    std::vector<double> criterion_variables(num_solution);

    for (uint i = 0; i < num_solution; ++i)
        criterion_variables[i] = P0_solutions[i][This_::criterion_variable_index_];

    const auto num_vnode = This_::vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_solution;
    vnode_index_to_allowable_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : This_::vnode_index_to_share_cell_indexes_) {
        const auto num_share_cell = share_cell_indexes.size();
        std::vector<double> target_solutions;
        target_solutions.reserve(num_share_cell);

        for (const auto cell_index : share_cell_indexes)
            target_solutions.push_back(criterion_variables[cell_index]);

        const auto min_solution = *std::min_element(target_solutions.begin(), target_solutions.end());
        const auto max_solution = *std::max_element(target_solutions.begin(), target_solutions.end());

        vnode_index_to_allowable_min_max_solution.emplace(vnode_index, std::make_pair(min_solution, max_solution));
    }

    return vnode_index_to_allowable_min_max_solution;
}

template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
bool hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::is_satisfy_P1_projected_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) {
    return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
}

template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
bool hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) {
    if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
        return true;
    else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
        return true;
    else
        return false;
}