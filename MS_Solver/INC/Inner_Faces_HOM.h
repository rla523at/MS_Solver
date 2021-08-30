#pragma once
#include "Grid_Builder.h"
#include "Reconstruction_Method_HOM.h"
#include "Numerical_Flux_Function.h"
#include "Solution_Scaling_Method.h"

//HOM�̸� �������� ����ϴ� variable & method
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
class Inner_Faces_HOM
{
private:
    static constexpr ushort space_dimension_    = Reconstruction_Method::space_dimension();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();
    static constexpr ushort num_equation_       = Numerical_Flux_Function::num_equation();
    
    using This_                 = Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;

protected:
    const Reconstruction_Method& reconstruction_method_;
    std::vector<std::pair<uint, uint>> oc_nc_index_pairs_;
    std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_qnodes_pairs_;
    std::vector<std::vector<Space_Vector_>> set_of_normals_;
    std::vector<std::pair<Dynamic_Matrix, Dynamic_Matrix>> oc_nc_side_basis_weight_pairs_;

public:
    Inner_Faces_HOM(Grid<space_dimension_>&& grid, const Reconstruction_Method& reconstruction_method);

    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const;

public:
    void initialize_scaling_method(void) const;
};


// template definition part
template<typename Reconstruction_Method, typename Numerical_Flux_Function>
Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>::Inner_Faces_HOM(Grid<space_dimension_>&& grid, const Reconstruction_Method& reconstruction_method)
    : reconstruction_method_(reconstruction_method){
    SET_TIME_POINT;

    this->oc_nc_index_pairs_ = std::move(grid.connectivity.inner_face_oc_nc_index_pairs);

    const auto& cell_elements = grid.elements.cell_elements;
    const auto& inner_face_elements = grid.elements.inner_face_elements;

    const auto num_inner_face = inner_face_elements.size();
    this->oc_nc_side_basis_qnodes_pairs_.reserve(num_inner_face);
    this->set_of_normals_.reserve(num_inner_face);
    this->oc_nc_side_basis_weight_pairs_.reserve(num_inner_face);

    constexpr auto integrand_order = 2 * Reconstruction_Method::solution_order() + 1;

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto& inner_face_element = inner_face_elements[i];
        const auto& inner_face_geometry = inner_face_element.geometry_;

        const auto& quadrature_rule = inner_face_geometry.get_quadrature_rule(integrand_order);
        const auto& qnodes = quadrature_rule.points;
        const auto& qweights = quadrature_rule.weights;

        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& oc_element = cell_elements[oc_index];

        auto oc_side_basis_qnode = this->reconstruction_method_.calculate_basis_nodes(oc_index, qnodes);
        auto nc_side_basis_qnode = this->reconstruction_method_.calculate_basis_nodes(nc_index, qnodes);

        this->oc_nc_side_basis_qnodes_pairs_.push_back({ std::move(oc_side_basis_qnode), std::move(nc_side_basis_qnode) });

        const auto num_qnode = qnodes.size();
        std::vector<Space_Vector_> normals(num_qnode);
        Dynamic_Matrix oc_side_basis_weight(num_qnode, This_::num_basis_);
        Dynamic_Matrix nc_side_basis_weight(num_qnode, This_::num_basis_);

        for (ushort q = 0; q < num_qnode; ++q) {
            normals[q] = inner_face_element.normalized_normal_vector(oc_element, qnodes[q]);
            oc_side_basis_weight.change_row(q, this->reconstruction_method_.calculate_basis_node(oc_index, qnodes[q]) * qweights[q]);          
            nc_side_basis_weight.change_row(q, this->reconstruction_method_.calculate_basis_node(nc_index, qnodes[q]) * qweights[q]); 
        }

        this->set_of_normals_.push_back(std::move(normals));
        this->oc_nc_side_basis_weight_pairs_.push_back({ std::move(oc_side_basis_weight), std::move(nc_side_basis_weight) });
    }

    Log::content_ << std::left << std::setw(50) << "@ Inner faces HOM precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}


template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const {
    const auto num_inner_face = this->set_of_normals_.size();

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];

        const auto& oc_solution_coefficient = solution_coefficients[oc_index];
        const auto& nc_solution_coefficient = solution_coefficients[nc_index];

        const auto& [oc_side_basis_qnode, nc_side_basis_qnode] = this->oc_nc_side_basis_qnodes_pairs_[i];
        const auto oc_side_cvariables = oc_solution_coefficient * oc_side_basis_qnode;
        const auto nc_side_cvariables = nc_solution_coefficient * nc_side_basis_qnode;

        const auto [num_equation, num_qnode] = oc_side_cvariables.size();
        const auto& normals = this->set_of_normals_[i];

        Dynamic_Matrix numerical_flux_quadrature(This_::num_equation_, num_qnode);
        for (ushort q = 0; q < num_qnode; ++q) {
            const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
            const auto nc_side_cvariable = nc_side_cvariables.column<This_::num_equation_>(q);

            numerical_flux_quadrature.change_column(q, Numerical_Flux_Function::calculate(oc_side_cvariable, nc_side_cvariable, normals[q]));
        }

        Residual_ owner_side_delta_rhs, neighbor_side_delta_rhs;
        const auto& [oc_side_basis_weight, nc_side_basis_weight] = this->oc_nc_side_basis_weight_pairs_[i];
        ms::gemm(numerical_flux_quadrature, oc_side_basis_weight, owner_side_delta_rhs);
        ms::gemm(numerical_flux_quadrature, nc_side_basis_weight, neighbor_side_delta_rhs);

        RHS[oc_index] -= owner_side_delta_rhs;
        RHS[nc_index] += neighbor_side_delta_rhs;
    }
}

template<typename Reconstruction_Method, typename Numerical_Flux_Function>
void Inner_Faces_HOM<Reconstruction_Method, Numerical_Flux_Function>::initialize_scaling_method(void) const {
    const auto num_inner_face = this->oc_nc_index_pairs_.size();

    for (uint i = 0; i < num_inner_face; ++i) {
        const auto [oc_index, nc_index] = this->oc_nc_index_pairs_[i];
        const auto& [oc_side_basis_qnodes, nc_side_basis_qnodes] = this->oc_nc_side_basis_qnodes_pairs_[i];

        Solution_Scaler::record_face_basis_qnodes(oc_index, oc_side_basis_qnodes);
        Solution_Scaler::record_face_basis_qnodes(nc_index, nc_side_basis_qnodes);
    }
}