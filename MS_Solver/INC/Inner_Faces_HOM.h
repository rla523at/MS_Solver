#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method.h"

//HOM이면 공통으로 사용하는 variable & method
template <typename Reconstruction_Method>
class Inner_Faces_HOM
{
private:
    static constexpr ushort space_dimension_ = Reconstruction_Method::space_dimension();
    static constexpr ushort num_basis_ = Reconstruction_Method::num_basis();
    
    using This_             = Inner_Faces_HOM<Reconstruction_Method>;
    using Space_Vector_     = Euclidean_Vector<space_dimension_>;    

protected:
    inline static std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    inline static std::vector<Dynamic_Matrix> set_of_oc_side_basis_qnodes_;
    inline static std::vector < std::vector<Space_Vector_>> set_of_normals_;
    inline static std::vector<Dynamic_Matrix> basis_weights_;

private:
    Inner_Faces_HOM(void) = delete;

public:
    static void initialize(Grid<space_dimension_>&& grid);
    //static void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) 

};


// template definition part
template<typename Reconstruction_Method>
void Inner_Faces_HOM<Reconstruction_Method>::initialize(Grid<space_dimension_>&& grid) {
    SET_TIME_POINT;

    This_::oc_nc_index_pairs_ = std::move(grid.connectivity.inner_face_oc_nc_index_pairs);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& inner_face_elements = grid.elements.inner_face_elements;

    const auto num_inner_face = inner_face_elements.size();
    This_::set_of_oc_side_basis_qnodes_.reserve(num_inner_face);
    This_::set_of_normals_.reserve(num_inner_face);
    This_::basis_weights_.reserve(num_inner_face);

    constexpr auto integrand_order = 2 * Reconstruction_Method::solution_order() + 1;

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto& inner_face_element = inner_face_elements[i];
        const auto& inner_face_geometry = inner_face_element.geometry_;

        const auto [oc_index, nc_index] = This_::oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        const auto& quadrature_rule = inner_face_geometry.get_quadrature_rule(integrand_order);
        const auto& qnodes = quadrature_rule.points;
        const auto& qweights = quadrature_rule.weights;
        const auto num_qnode = qnodes.size();

        This_::set_of_oc_side_basis_qnodes_.push_back(Reconstruction_Method::calculate_basis_nodes(oc_index, qnodes));

        std::vector<Space_Vector_> normals(num_qnode);
        Dynamic_Matrix basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            normals[q] = inner_face_element.normalized_normal_vector(oc_element, qnodes[q]);
            basis_weight.change_row(q, Reconstruction_Method::calculate_basis_node(oc_index, qnodes[q]) * qweights[q]);
        }

        This_::set_of_normals_.push_back(std::move(normals));
        This_::basis_weights_.push_back(std::move(basis_weight));
    }

    Log::content_ << std::left << std::setw(50) << "@ Inner faces FVM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template<typename Reconstruction_Method>
void Inner_Faces_HOM<Reconstruction_Method>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) {
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