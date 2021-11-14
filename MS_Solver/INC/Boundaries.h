#pragma once
#include "Boundary_Flux_Function.h"
#include "Discrete_Solution.h"
#include "Residual.h"


class Boundaries_DG
{
public:
    Boundaries_DG(const Configuration& configuration, const Grid& grid, const Discrete_Solution_DG& discrete_solution)
    {
        SET_TIME_POINT;

        this->oc_indexes_ = grid.boundary_owner_cell_indexes();
        this->boundary_flux_functions_ = grid.boundary_flux_functions<Numerical_Flux_Function>();

        constexpr auto integrand_degree = 2 * Reconstruction_Method::solution_order() + 1;
        auto boundary_quadrature_rules = grid.boundary_quadrature_rules(integrand_degree);

        const auto num_boundary = this->oc_indexes_.size();
        this->set_of_oc_side_basis_qnodes_.reserve(num_boundary);
        this->set_of_oc_side_qweights_basis_.reserve(num_boundary);

        std::vector<std::vector<Euclidean_Vector<space_dimension_>>> set_of_qnodes(num_boundary);

        for (uint i = 0; i < num_boundary; ++i) {
            auto& quadrature_rule = boundary_quadrature_rules[i];
            auto& qnodes = quadrature_rule.points;
            const auto& qweights = quadrature_rule.weights;
            const auto num_qnode = qnodes.size();

            auto basis_qnodes = reconstruction_method.basis_nodes(this->oc_indexes_[i], qnodes);

            Matrix qweight_basis(num_qnode, This_::num_basis_);
            for (ushort q = 0; q < num_qnode; ++q)
                qweight_basis.change_row(q, reconstruction_method.calculate_basis_node(this->oc_indexes_[i], qnodes[q]) * qweights[q]);

            this->set_of_oc_side_basis_qnodes_.push_back(std::move(basis_qnodes));
            this->set_of_oc_side_qweights_basis_.push_back(std::move(qweight_basis));
            set_of_qnodes[i] = std::move(qnodes);
        }

        this->set_of_normals_ = grid.boundary_set_of_normals(this->oc_indexes_, set_of_qnodes);

        Log::content_ << std::left << std::setw(50) << "@ Boundaries HOM base precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
        Log::print();
    }

public:
    void calculate_RHS(Residual& residual, const Discrete_Solution_DG& discrete_soltuion) const
    { 
        for (int i = 0; i < this->num_boundaries_; ++i) 
        {
            const auto oc_index = this->oc_indexes_[i];
            const auto solution_at_QPs = discrete_soltuion.calculate_solution_at_points(oc_index, this->set_of_oc_side_basis_QPs_m_[i]);
            const auto num_QPs = solution_at_QPs.size();

            const auto& boundary_flux_function = *this->boundary_flux_functions_[i];
            const auto& normals = this->set_of_normals_[i];
            
            Matrix boundary_flux_quadrature(this->num_equations_, num_QPs);

            for (ushort q = 0; q < num_QPs; ++q) 
            {
                boundary_flux_quadrature.change_column(q, boundary_flux_function.calculate(solution_at_QPs[q], normals[q]));
            }

            const auto delta_rhs = boundary_flux_quadrature * this->set_of_oc_side_QWs_basis_m_[i]; // RHS[oc_index] -= owner_side_delta_rhs; 만들때 -붙혀서 만들기!
            residual.update_rhs(oc_index, delta_rhs);            
        }
    }

private:
    uint num_boundaries_;
    ushort num_equations_;
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