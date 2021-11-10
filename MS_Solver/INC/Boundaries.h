#pragma once
#include "Boundary_Flux_Function.h"
#include "Discrete_Solution.h"

class Boundaries_DG
{
public:
    void calculate_RHS(const double* RHS, const Discrete_Solution_DG& discrete_soltuion) const
    { 
        for (int i = 0; i < this->num_boundaries_; ++i) {
            const auto oc_index = this->oc_indexes_[i];
            const auto solution_at_QPs = discrete_soltuion.calculate_solution_at_points(oc_index, this->set_of_oc_side_basis_QPs_m_[i]);

            const auto oc_side_cvariables = this->set_of_oc_coefficient_m_[i] * this->set_of_oc_side_basis_QPs_m_[i];

            const auto [num_equation, num_qnode] = oc_side_cvariables.size();

            const auto& boundary_flux_function = *this->boundary_flux_functions_[i];
            const auto& normals = this->set_of_normals_[i];

            Matrix boundary_flux_quadrature(This_::num_equation_, num_qnode);
            for (ushort q = 0; q < num_qnode; ++q) {
                const auto oc_side_cvariable = oc_side_cvariables.column<This_::num_equation_>(q);
                boundary_flux_quadrature.change_column(q, boundary_flux_function.calculate(oc_side_cvariable, normals[q]));
            }

            Residual_ owner_side_delta_rhs;
            ms::gemm(boundary_flux_quadrature, this->set_of_oc_side_qweights_basis_[i], owner_side_delta_rhs);

            RHS[oc_index] -= owner_side_delta_rhs;
        }
    }

private:
    uint num_boundaries_;
    std::vector<std::unique_ptr<Boundary_Flux_Function>> boundary_flux_functions_;
    std::vector<uint> oc_indexes_;
    std::vector<Matrix> set_of_oc_side_basis_QPs_m_;
    std::vector<Matrix> set_of_oc_side_QWs_basis_m_;
    std::vector<std::vector<Euclidean_Vector>> set_of_normals_;
};



//#include "Boundaries_FVM.h"
//#include "Boundaries_HOM.h"
//#include "Spatial_Discrete_Method.h"
//
//template <typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Boundaries;
//
//
//template <typename Numerical_Flux_Function>
//class Boundaries<FVM, Constant_Reconstruction, Numerical_Flux_Function> : public Boundaries_FVM_Constant<Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Boundaries(const Grid<space_dimension_>& grid, const Constant_Reconstruction& reconstruction_method)
//        : Boundaries_FVM_Constant<Numerical_Flux_Function>(grid) {};
//};
//
//
//template <typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Boundaries<FVM, Reconstruction_Method, Numerical_Flux_Function> : public Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Boundaries(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
//        : Boundaries_FVM_Linear<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
//};
//
//
//template <typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Boundaries<HOM, Reconstruction_Method, Numerical_Flux_Function> : public Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>
//{
//private:
//    static constexpr size_t space_dimension_ = Numerical_Flux_Function::space_dimension();
//
//public:
//    Boundaries(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method)
//        : Boundaries_HOM<Reconstruction_Method, Numerical_Flux_Function>(grid, reconstruction_method) {};
//};