#pragma once
#include "Boundary_Flux_Function.h"
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"


//HOM이면 공통으로 사용하는 variable
template <typename Governing_Equation, typename Reconstruction_Method>
class Boundaries_HOM
{
private:
    static constexpr ushort space_dimension_    = Governing_Equation::space_dimension();
    static constexpr ushort num_equation_       = Governing_Equation::num_equation();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();

    using This_                 = Boundaries_HOM<Governing_Equation, Reconstruction_Method>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;

private:
    Boundaries_HOM(void) = delete;

protected:
    inline static std::vector<std::unique_ptr<Boundary_Flux_Function<Governing_Equation>>> boundary_flux_functions_;

    inline static std::vector<uint> oc_indexes_;
    inline static std::vector<Dynamic_Matrix> set_of_oc_side_basis_qnodes_;
    inline static std::vector<std::vector<Space_Vector_>> set_of_normals_;

    inline static std::vector<Dynamic_Matrix> basis_weights_;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients);
};


//template definition
template <typename Governing_Equation, typename Reconstruction_Method>
void Boundaries_HOM<Governing_Equation, Reconstruction_Method>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    This_::oc_indexes_ = std::move(grid.connectivity.boundary_oc_indexes);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& boundary_elements = grid.elements.boundary_elements;

    const auto num_boundary = boundary_elements.size();
    This_::boundary_flux_functions_.reserve(num_boundary);
    This_::set_of_oc_side_basis_qnodes_.reserve(num_boundary);
    This_::set_of_normals_.reserve(num_boundary);
    This_::basis_weights_.reserve(num_boundary);

    constexpr auto integrand_order = 2 * Reconstruction_Method::solution_order() + 1;

    for (uint i = 0; i < num_boundary; ++i) {
        const auto& boundary_element = boundary_elements[i];
        const auto& boundary_geometry = boundary_element.geometry_;

        This_::boundary_flux_functions_.push_back(Boundary_Flux_Function_Factory<Governing_Equation>::make(boundary_element.type()));

        const auto oc_index = This_::oc_indexes_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto& quadrature_rule = boundary_geometry.get_quadrature_rule(integrand_order);
        const auto& qnodes = quadrature_rule.points;
        const auto& qweights = quadrature_rule.weights;
        const auto num_qnode = qnodes.size();

        This_::set_of_oc_side_basis_qnodes_.push_back(Reconstruction_Method::calculate_basis_nodes(oc_index, qnodes));

        std::vector<Space_Vector_> normals(num_qnode);
        Dynamic_Matrix basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            normals[q] = boundary_element.normalized_normal_vector(oc_element, qnodes[q]);
            basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(oc_index, qnodes[q]) * qweights[q]);
        }

        This_::set_of_normals_.push_back(std::move(normals));
        This_::basis_weights_.push_back(std::move(basis_weight));
    }

    Log::content_ << std::left << std::setw(50) << "@ Boundaries FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template <typename Governing_Equation, typename Reconstruction_Method>
void Boundaries_HOM<Governing_Equation, Reconstruction_Method>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) {
    const auto num_boundary = This_::boundary_flux_functions_.size();

    for (uint i = 0; i < num_boundary; ++i) {
        const auto oc_index = This_::oc_indexes_[i];

        const auto& solution_coefficient = solution_coefficients[oc_index];
        const auto oc_side_cvariables = solution_coefficient * This_::set_of_oc_side_basis_qnodes_[i];

        const auto [num_equation, num_qnode] = oc_side_cvariables.size();

        const auto& boundary_flux_function = *This_::boundary_flux_functions_[i];
        const auto& normals = This_::set_of_normals_[i];

        Dynamic_Matrix numerical_flux_quadrature(This_::num_equation_, num_qnode);
        for (ushort q = 0; q < num_qnode; ++q) {
            const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
            numerical_flux_quadrature.change_column(q, boundary_flux_function.calculate(oc_side_cvariable, normals[q]));
        }

        ms::gemm(numerical_flux_quadrature, basis_weights_[i], RHS[oc_index]);
    }
}