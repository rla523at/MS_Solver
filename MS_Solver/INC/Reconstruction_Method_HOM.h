#pragma once
#include "Governing_Equation.h"
#include "Reconstruction_Method_FVM.h"
#include "Polynomial.h"


class HOM_Reconstruction : public RM {};

template <ushort space_dimension_, ushort solution_order_>
class Polynomial_Reconstruction : public HOM_Reconstruction
{
protected:
    using Space_Vector_ = Euclidean_Vector<space_dimension_>;

    static constexpr ushort num_basis_ = ms::combination_with_repetition(1 + space_dimension_, solution_order_);

protected:
    std::vector<Vector_Function<Polynomial<space_dimension_>, num_basis_>> set_of_basis_functions_;

public:
    Polynomial_Reconstruction(const Grid<space_dimension_>& grid);

    auto calculate_set_of_transposed_gradient_basis(void) const;
    auto calculate_basis_node(const uint cell_index, const Space_Vector_& node) const;
    Dynamic_Matrix calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) const;
    double calculate_P0_basis_value(const uint cell_index, const Space_Vector_& node) const;

    //std::vector<Dynamic_Matrix> calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes) const;
    //std::vector<double> calculate_P0_basis_values(const std::vector<Space_Vector_>& nodes) const;

public:
    static constexpr ushort num_basis(void) { return num_basis_; };
    static constexpr ushort solution_order(void) { return solution_order_; };
    static constexpr ushort space_dimension(void) { return space_dimension_; };
    static std::string name(void) { return "P" + std::to_string(solution_order_) + "_Polynomial_Reconstruction"; };
};


class P1_Projected_MLP_Condition
{
public:
    static bool is_satisfy(const double P1_projected_value, const double allowable_min, const double allowable_max) {
        return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;;
    }
};

class MLP_Smooth_Extrema_Detector
{
public:
    static bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) {
        if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
            return true;
        else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
            return true;
        else
            return false;
    }
};

template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
class hMLP_Reconstruction : public Polynomial_Reconstruction<space_dimension_, solution_order_>
{
private:
    static_require(1 <= solution_order_, "solution order should be greater than 1 for hMLP reconstruction");
    
    static constexpr ushort num_basis_ = ms::combination_with_repetition(1 + space_dimension_, solution_order_);
    static constexpr ushort criterion_variable_index_ = 0; //rho

    using This_                     = hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>;
    using Solution_Coefficients_    = Matrix<num_equation, num_basis_>;

private:
    std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_indexes_;
    std::vector<std::vector<uint>> set_of_vnode_indexes_;
    std::array<ushort, solution_order_ + 1> num_Pn_projection_basis_;
    std::vector<Dynamic_Matrix> set_of_basis_vnodes_;
    std::vector<Dynamic_Matrix> set_of_P1_projected_basis_vnodes_;
    std::vector<Dynamic_Matrix> set_of_P1_mode_basis_vnodes_;
    std::vector<double> P0_basis_values_;    
    
public:
    hMLP_Reconstruction(Grid<space_dimension_>&& grid);

    void reconstruct(std::vector<Solution_Coefficients_>& solutions);
    
private:
    auto calculate_vertex_node_index_to_allowable_min_max_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions) const;
    auto Pn_projection_matrix(const ushort Pn) const;

public:
    static std::string name(void) { return "P" + std::to_string(solution_order_) + "_hMLP_Reconstruction"; };
};


//template definition part
template <ushort space_dimension_, ushort solution_order_>
Polynomial_Reconstruction<space_dimension_, solution_order_>::Polynomial_Reconstruction(const Grid<space_dimension_>& grid) {
    SET_TIME_POINT;

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->set_of_basis_functions_.reserve(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        this->set_of_basis_functions_.push_back(cell_geometry.orthonormal_basis_functions<solution_order_>());
    }

    Log::content_ << std::left << std::setw(50) << "@ Polynomial reconstruction precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_transposed_gradient_basis(void) const {
    const auto num_cell = this->set_of_basis_functions_.size();
    std::vector<Matrix_Function<Polynomial<space_dimension_>, space_dimension_, this->num_basis_>> set_of_transposed_gradient_basis(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        auto& transposed_gradient_basis = set_of_transposed_gradient_basis[i];
        const auto& basis_function = this->set_of_basis_functions_[i];

        for (ushort j = 0; j < this->num_basis_; ++j)
            transposed_gradient_basis.change_column(j, basis_function[j].gradient());
    }

    return set_of_transposed_gradient_basis;
}

template <ushort space_dimension_, ushort solution_order_>
auto Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_node(const uint cell_index, const Space_Vector_& node) const {
    const auto& basis_functions = this->set_of_basis_functions_[cell_index];
    return basis_functions(node);
}


template <ushort space_dimension_, ushort solution_order_>
Dynamic_Matrix Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_basis_nodes(const uint cell_index, const std::vector<Space_Vector_>& nodes) const {
    const auto& basis_functions = this->set_of_basis_functions_[cell_index];

    const auto num_node = nodes.size();
    Dynamic_Matrix basis_nodes(this->num_basis_, num_node);

    for (ushort j = 0; j < num_node; ++j)
        basis_nodes.change_column(j, basis_functions(nodes[j]));

    return basis_nodes;
}

template <ushort space_dimension_, ushort solution_order_>
double Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_value(const uint cell_index, const Space_Vector_& node) const {
    const auto& basis_function = this->set_of_basis_functions_[cell_index];
    const auto& P0_basis_function = basis_function[0];

    return P0_basis_function(node);
}

//template <ushort space_dimension_, ushort solution_order_>
//std::vector<Dynamic_Matrix> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_set_of_basis_nodes(const std::vector<std::vector<Space_Vector_>>& set_of_nodes) const {
//    const auto num_cell = this->set_of_basis_functions_.size();
//
//    dynamic_require(set_of_nodes.size() == num_cell, "set size should be same with num cell");
//
//    std::vector<Dynamic_Matrix> set_of_basis_nodes;
//    set_of_basis_nodes.reserve(num_cell);
//
//    for (uint i = 0; i < num_cell; ++i)
//        set_of_basis_nodes.push_back(this->calculate_basis_nodes(i, set_of_nodes[i]));
//
//    return set_of_basis_nodes;
//}

//template <ushort space_dimension_, ushort solution_order_>
//std::vector<double> Polynomial_Reconstruction<space_dimension_, solution_order_>::calculate_P0_basis_values(const std::vector<Space_Vector_>& nodes) const {
//    const auto num_nodes = nodes.size();
//    std::vector<double> P0_basis_values(num_nodes);
//
//    for (uint i = 0; i < num_nodes; ++i) {
//        const auto& basis_function = this->set_of_basis_functions_[i];
//        const auto& P0_basis_function = basis_function[0];
//        const auto& cell_center = nodes[i];
//
//        P0_basis_values[i] = P0_basis_function(cell_center);
//    }
//
//    return P0_basis_values;
//}


template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::hMLP_Reconstruction(Grid<space_dimension_>&& grid) 
    : Polynomial_Reconstruction<space_dimension_, solution_order_>(grid){
    SET_TIME_POINT;
    
    for (ushort i = 0; i <= solution_order_; ++i) 
        this->num_Pn_projection_basis_[i] = ms::combination_with_repetition(1 + space_dimension_, i);

    this->vnode_index_to_share_cell_indexes_ = std::move(grid.connectivity.vnode_index_to_share_cell_indexes);

    const auto& cell_elements = grid.elements.cell_elements;

    const auto num_cell = cell_elements.size();
    this->set_of_vnode_indexes_.reserve(num_cell);
    this->set_of_basis_vnodes_.reserve(num_cell);
    this->set_of_P1_projected_basis_vnodes_.reserve(num_cell);
    this->set_of_P1_mode_basis_vnodes_.reserve(num_cell);
    this->P0_basis_values_.reserve(num_cell);

    constexpr ushort P0_basis_row_index = 0;

    for (uint i = 0; i < num_cell; ++i) {
        const auto& cell_element = cell_elements[i];
        const auto& cell_geometry = cell_element.geometry_;

        this->set_of_vnode_indexes_.push_back(cell_element.vertex_node_indexes());

        const auto vnodes = cell_geometry.vertex_nodes();
        const auto num_vnode = vnodes.size();

        auto basis_vnodes = this->calculate_basis_nodes(i, vnodes);

        const auto P1_projection_matrix = this->Pn_projection_matrix(1);
        auto P1_projected_basis_vnodes = P1_projection_matrix * basis_vnodes;

        Dynamic_Matrix P1_mode_basis_vnodes = P1_projected_basis_vnodes;
        P1_mode_basis_vnodes.change_row(P0_basis_row_index, Dynamic_Euclidean_Vector(num_vnode));

        const auto P0_basis_value = P1_projected_basis_vnodes.at(P0_basis_row_index, 0);

        this->set_of_basis_vnodes_.push_back(std::move(basis_vnodes));
        this->set_of_P1_projected_basis_vnodes_.push_back(std::move(P1_projected_basis_vnodes));
        this->set_of_P1_mode_basis_vnodes_.push_back(std::move(P1_mode_basis_vnodes));
        this->P0_basis_values_.push_back(P0_basis_value);
    }

    Log::content_ << std::left << std::setw(50) << "@ hMLP precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <ushort num_equation, ushort space_dimension_, ushort solution_order_>
void hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::reconstruct(std::vector<Solution_Coefficients_>& solution_coefficients) {
    const auto num_cell = solution_coefficients.size();
    std::vector<Euclidean_Vector<num_equation>> P0_solutions(num_cell);

    for (uint i = 0; i < num_cell; ++i) {
        const auto P0_coefficient = solution_coefficients[i].column(0);
        P0_solutions[i] = P0_coefficient * this->P0_basis_values_[i];
    }

    const auto vnode_index_to_allowable_min_max_solution = this->calculate_vertex_node_index_to_allowable_min_max_solution(P0_solutions);

    for (uint i = 0; i < num_cell; ++i) {
        auto temporal_solution_order = solution_order_;

        auto& solution_coefficient = solution_coefficients[i];
        const auto& basis_vnodes = this->set_of_basis_vnodes_[i];
        const auto& P1_projected_basis_vnodes = this->set_of_P1_projected_basis_vnodes_[i];

        auto solution_vnodes = solution_coefficient * basis_vnodes;
        auto P1_projected_solution_vnodes = solution_coefficient * P1_projected_basis_vnodes;

        const auto& vnode_indexes = this->set_of_vnode_indexes_[i];
        const auto num_vnode = vnode_indexes.size();

        bool end_limiting = false;

        for (ushort j = 0; j < num_vnode; ++j) {
            const auto vnode_index = vnode_indexes[j];
            const auto [allowable_min, allowable_max] = vnode_index_to_allowable_min_max_solution.at(vnode_index);

            while (true) {
                const auto P1_projected_solution = P1_projected_solution_vnodes.column<num_equation>(j);
                const auto P1_projected_solution_criterion_variable = P1_projected_solution[This_::criterion_variable_index_];

                if (P1_Projected_MLP_Condition::is_satisfy(P1_projected_solution_criterion_variable, allowable_min, allowable_max))
                    break;
                else {
                    const auto solution = solution_vnodes.column<num_equation>(j);

                    const auto solution_criterion_variable = solution[This_::criterion_variable_index_];
                    const auto higher_mode_criterion_variable = solution_criterion_variable - P1_projected_solution_criterion_variable;
                    const auto P1_mode_criterion_variable = P1_projected_solution_criterion_variable - P0_solutions[i][This_::criterion_variable_index_];

                    if (MLP_Smooth_Extrema_Detector::is_smooth_extrema(solution_criterion_variable, higher_mode_criterion_variable, P1_mode_criterion_variable, allowable_min, allowable_max))
                        break;
                    else {
                        if (temporal_solution_order == 1) {
                            const auto P0_mode_criterion_variable = P0_solutions[i][This_::criterion_variable_index_];

                            const auto limiting_value = MLP_u1_Limiting_Strategy::limit(P1_mode_criterion_variable, P0_mode_criterion_variable, allowable_min, allowable_max);
                            
                            std::array<double, This_::num_basis_> limiting_values;
                            limiting_values.fill(limiting_value);
                            limiting_values[0] = 1.0; // keep P0 coefficient
                            
                            const Matrix limiting_matrix = limiting_values;
                            solution_coefficient *= limiting_matrix;

                            end_limiting = true;
                            break;
                        }
                        else {
                            //limiting highest mode
                            const auto Pnm1_projection_matrix = this->Pn_projection_matrix(--temporal_solution_order);
                            solution_coefficient *= Pnm1_projection_matrix;

                            //re-calculate limited vnode_solution                            
                            solution_vnodes = solution_coefficient * basis_vnodes;
                            P1_projected_solution_vnodes = solution_coefficient * P1_projected_basis_vnodes;
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
auto hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::calculate_vertex_node_index_to_allowable_min_max_solution(const std::vector<Euclidean_Vector<num_equation>>& P0_solutions) const {
    const auto num_solution = P0_solutions.size();
    std::vector<double> criterion_variables(num_solution);

    for (uint i = 0; i < num_solution; ++i)
        criterion_variables[i] = P0_solutions[i][this->criterion_variable_index_];

    const auto num_vnode = this->vnode_index_to_share_cell_indexes_.size();

    std::unordered_map<uint, std::pair<double, double>> vnode_index_to_allowable_min_max_solution;
    vnode_index_to_allowable_min_max_solution.reserve(num_vnode);

    for (const auto& [vnode_index, share_cell_indexes] : this->vnode_index_to_share_cell_indexes_) {
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
auto hMLP_Reconstruction<num_equation, space_dimension_, solution_order_>::Pn_projection_matrix(const ushort Pn) const {
    dynamic_require(Pn <= solution_order_, "Projection order should be less then solution order");

    const auto num_Pn_projection_basis = this->num_Pn_projection_basis_[Pn];
    
    std::array<double, this->num_basis_> limiting_value = { 0 };
    for (ushort j = 0; j < num_Pn_projection_basis; ++j)
        limiting_value[j] = 1.0; //preserve Pn_projection_basis

    Matrix Pn_projection_matrix = limiting_value;
    return Pn_projection_matrix;
 }