#pragma once
#include "Boundary_Flux_Function.h"
#include "Grid.h"
#include "Reconstruction_Method_HOM.h"
#include "Solution_Scaler.h"


//HOM이면 공통으로 사용하는 variable
template <typename Reconstruction_Method, typename Numerical_Flux_Function>
class Boundaries_HOM
{
private:
    static constexpr ushort space_dimension_    = Numerical_Flux_Function::space_dimension();
    static constexpr ushort num_equation_       = Numerical_Flux_Function::num_equation();
    static constexpr ushort num_basis_          = Reconstruction_Method::num_basis();

    using This_                 = Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>;
    using Space_Vector_         = Euclidean_Vector<space_dimension_>;
    using Residual_             = Matrix<num_equation_, num_basis_>;
    using Solution_Coefficient_ = Matrix<num_equation_, num_basis_>;

protected:
    const Reconstruction_Method& reconstruction_method_;
    std::vector<std::unique_ptr<Boundary_Flux_Function<Numerical_Flux_Function>>> boundary_flux_functions_;
    std::vector<uint> oc_indexes_;
    std::vector<Dynamic_Matrix> set_of_oc_side_basis_qnodes_;
    std::vector<std::vector<Space_Vector_>> set_of_normals_;
    std::vector<Dynamic_Matrix> set_of_oc_side_basis_weights_;

public:
    Boundaries_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method);

public:
    void calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const;

public:
    void initialize_scaling_method(void) const;
};


//template definition
template <typename Reconstruction_Method, typename Numerical_Flux_Function>
Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::Boundaries_HOM(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
    : reconstruction_method_(reconstruction_method) {
    SET_TIME_POINT;

    this->oc_indexes_ = grid.boundary_owner_cell_indexes();
    this->boundary_flux_functions_ = grid.boundary_flux_functions<Numerical_Flux_Function>();

    constexpr auto integrand_degree = 2 * Reconstruction_Method::solution_order() + 1;
    auto boundary_quadrature_rules = grid.boundary_quadrature_rules(integrand_degree);

    const auto num_boundary = this->oc_indexes_.size();
    this->set_of_oc_side_basis_qnodes_.reserve(num_boundary);
    this->set_of_oc_side_basis_weights_.reserve(num_boundary);

    std::vector<std::vector<Euclidean_Vector<space_dimension_>>> set_of_qnodes;
    set_of_qnodes.reserve(num_boundary);

    for (uint i = 0; i < num_boundary; ++i) {
        auto& quadrature_rule = boundary_quadrature_rules[i];
        auto& qnodes = quadrature_rule.points;
        const auto& qweights = quadrature_rule.weights;
        const auto num_qnode = qnodes.size();

        Dynamic_Matrix basis_weights(num_qnode, This_::num_basis_);
        for (ushort q = 0; q < num_qnode; ++q) 
            basis_weights.change_row(q, reconstruction_method.calculate_basis_node(this->oc_indexes_[i], qnodes[q]) * qweights[q]);
        
        this->set_of_oc_side_basis_qnodes_.push_back(reconstruction_method.calculate_basis_nodes(this->oc_indexes_[i], qnodes));
        this->set_of_oc_side_basis_weights_.push_back(std::move(basis_weights));
        set_of_qnodes.push_back(std::move(qnodes));
    }

    this->set_of_normals_ = grid.boundary_set_of_normals(this->oc_indexes_, set_of_qnodes);

    Log::content_ << std::left << std::setw(50) << "@ Boundaries HOM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
    Log::print();
}

template <typename Reconstruction_Method, typename Numerical_Flux_Function>
void Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::calculate_RHS(std::vector<Residual_>& RHS, const std::vector<Solution_Coefficient_>& solution_coefficients) const {
    const auto num_boundary = this->boundary_flux_functions_.size();

    for (uint i = 0; i < num_boundary; ++i) {
        const auto oc_index = this->oc_indexes_[i];

        const auto& solution_coefficient = solution_coefficients[oc_index];
        const auto oc_side_cvariables = solution_coefficient * this->set_of_oc_side_basis_qnodes_[i];

        const auto [num_equation, num_qnode] = oc_side_cvariables.size();

        const auto& boundary_flux_function = *this->boundary_flux_functions_[i];
        const auto& normals = this->set_of_normals_[i];

        Dynamic_Matrix boundary_flux_quadrature(This_::num_equation_, num_qnode);
        for (ushort q = 0; q < num_qnode; ++q) {
            const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
            boundary_flux_quadrature.change_column(q, boundary_flux_function.calculate(oc_side_cvariable, normals[q]));
        }

        Residual_ owner_side_delta_rhs;
        ms::gemm(boundary_flux_quadrature,this->set_of_oc_side_basis_weights_[i], owner_side_delta_rhs);

        RHS[oc_index] -= owner_side_delta_rhs;
    }
}


template <typename Reconstruction_Method, typename Numerical_Flux_Function>
void Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>::initialize_scaling_method(void) const {
    const auto num_boundary = this->oc_indexes_.size();
    for (uint i = 0; i < num_boundary; ++i) 
        Solution_Scaler<This_::space_dimension_>::record_face_basis_qnodes(this->oc_indexes_[i], this->set_of_oc_side_basis_qnodes_[i]);    
}
