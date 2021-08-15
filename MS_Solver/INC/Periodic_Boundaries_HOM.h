#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method_HOM.h"
#include "Numerical_Flux_Function.h"

template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
    for (const auto& elem : vec)
        os << elem << "\n";
    return os;
}


//HOM이면 공통으로 사용하는 variable & method
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
class Periodic_Boundaries_HOM
{
private:
    static constexpr ushort space_dimension_    = Reconstruction_Method::space_dimension();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();
    static constexpr ushort num_equation_       = Numerical_Flux_Function::num_equation();

    using This_                 = Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;

protected:
    inline static std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    inline static std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_qnode_pairs_;
    inline static std::vector<std::vector<Space_Vector_>> set_of_normals_;
    inline static std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_weight_pairs_;

private:
    Periodic_Boundaries_HOM(void) = delete;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients);
};


// template definition part
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    This_::oc_nc_index_pairs_ = std::move(grid.connectivity.periodic_boundary_oc_nc_index_pairs);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& periodic_boundary_element_pairs = grid.elements.periodic_boundary_element_pairs;

    const auto num_periodic_pair = periodic_boundary_element_pairs.size();
    This_::oc_nc_side_basis_qnode_pairs_.reserve(num_periodic_pair);
    This_::set_of_normals_.reserve(num_periodic_pair);
    This_::oc_nc_side_basis_weight_pairs_.reserve(num_periodic_pair);

    constexpr auto integrand_order = 2 * Reconstruction_Method::solution_order() + 1;

    for (uint i = 0; i < num_periodic_pair; ++i) {
        const auto& [oc_side_element, nc_side_element] = periodic_boundary_element_pairs[i];
        const auto& oc_side_geometry = oc_side_element.geometry_;
        const auto& nc_side_geometry = nc_side_element.geometry_;

        const auto& oc_side_quadrature_rule = oc_side_geometry.get_quadrature_rule(integrand_order);
        const auto& oc_side_qnodes = oc_side_quadrature_rule.points;
        const auto& oc_side_qweights = oc_side_quadrature_rule.weights;

        const auto& nc_side_quadrature_rule = nc_side_geometry.get_quadrature_rule(integrand_order);
        const auto& nc_side_qnodes = nc_side_quadrature_rule.points;
        const auto& nc_side_qweights = nc_side_quadrature_rule.weights;

        //debug
        auto re_ordered_nc_side_qnodes = nc_side_qnodes;
        auto re_ordered_nc_side_qweights = nc_side_qweights;

        const auto oc_side_normal_vector = oc_side_geometry.normalized_normal_vector(oc_side_geometry.center_node());
        const auto nc_side_normal_vector = nc_side_geometry.normalized_normal_vector(nc_side_geometry.center_node());
        const auto inner_product = oc_side_normal_vector.inner_product(nc_side_normal_vector);

        if (inner_product < 0) {
            std::reverse(re_ordered_nc_side_qnodes.begin(), re_ordered_nc_side_qnodes.end());
            std::reverse(re_ordered_nc_side_qweights.begin(), re_ordered_nc_side_qweights.end());
        }   
        //debug

        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        auto oc_side_basis_qnode = Reconstruction_Method::calculate_basis_nodes(oc_index, oc_side_qnodes);
        auto nc_side_basis_qnode = Reconstruction_Method::calculate_basis_nodes(nc_index, re_ordered_nc_side_qnodes);

        const auto num_qnode = oc_side_qnodes.size();
        std::vector<Space_Vector_> normals(num_qnode);
        Dynamic_Matrix oc_side_basis_weight(num_qnode, This_::num_basis_);
        Dynamic_Matrix nc_side_basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            normals[q] = oc_side_element.normalized_normal_vector(oc_element, oc_side_qnodes[q]);
            oc_side_basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(oc_index, oc_side_qnodes[q]) * oc_side_qweights[q]);    
            nc_side_basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(nc_index, re_ordered_nc_side_qnodes[q]) * re_ordered_nc_side_qweights[q]);
        }
        

        This_::oc_nc_side_basis_qnode_pairs_.push_back({ std::move(oc_side_basis_qnode), std::move(nc_side_basis_qnode) });
        This_::set_of_normals_.push_back(std::move(normals));
        This_::oc_nc_side_basis_weight_pairs_.push_back({ std::move(oc_side_basis_weight), std::move(nc_side_basis_weight) });
    }

    Log::content_ << std::left << std::setw(50) << "@ Periodic Boundaries HOM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Periodic_Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) {
    const auto num_periodic_pairs = This_::oc_nc_index_pairs_.size();

    for (uint i = 0; i < num_periodic_pairs; ++i) {
        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];

        const auto& oc_solution_coefficient = solution_coefficients[oc_index];
        const auto& nc_solution_coefficient = solution_coefficients[nc_index];

        const auto& [oc_side_basis_qnode, nc_side_basis_qnode] = This_::oc_nc_side_basis_qnode_pairs_[i];
        const auto oc_side_cvariables = oc_solution_coefficient * oc_side_basis_qnode;
        const auto nc_side_cvariables = nc_solution_coefficient * nc_side_basis_qnode;

        const auto [num_equation, num_qnode] = oc_side_cvariables.size();
        const auto& normals = This_::set_of_normals_[i];

        Dynamic_Matrix numerical_flux_quadrature(This_::num_equation_, num_qnode);
        for (ushort q = 0; q < num_qnode; ++q) {
            const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
            const auto nc_side_cvariable = nc_side_cvariables.column<This_::num_equation_>(q);

            numerical_flux_quadrature.change_column(q, Numerical_Flux_Function::calculate(oc_side_cvariable, nc_side_cvariable, normals[q]));
        }

        This_::Residual_ owner_side_delta_rhs, neighbor_side_delta_rhs;
        const auto& [oc_side_basis_weight, nc_side_basis_weight] = This_::oc_nc_side_basis_weight_pairs_[i];
        ms::gemm(numerical_flux_quadrature, oc_side_basis_weight, owner_side_delta_rhs);
        ms::gemm(numerical_flux_quadrature, nc_side_basis_weight, neighbor_side_delta_rhs);

        RHS[oc_index] -= owner_side_delta_rhs;
        RHS[nc_index] += neighbor_side_delta_rhs;
    }
}