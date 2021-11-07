#pragma once
#include "Discrete_Solution.h"

class Cells
{
public:
    virtual void calculate_RHS(double* RHS, const Discrete_Solution& discretized_solution) abstract;
};


class Cells_HOM : public Cells
{
public:
    void calculate_RHS(double* RHS, const Discrete_Solution& discretized_solution) override
    {
        for (uint i = 0; i < this->num_cells_; ++i) {
            const auto solution_at_quadrature_nodes = discretized_solution.calculate_solution_at_quadrature_nodes(i);

            Matrix flux_quadrature_points(num_eq, This_::space_dimension_ * num_quadrature_node);

            for (size_t j = 0; j < num_quadrature_node; ++j) {
                const auto physical_flux = this->governing_equation_->calculate_physical_flux(solution_at_quadrature_nodes[j]);
                flux_quadrature_points.change_columns(j * This_::space_dimension_, physical_flux);
            }

            ms::gemm(flux_quadrature_points, This_::qweights_gradient_basis_[i], delta_rhs);

            RHS[i] += delta_rhs;
        }
    }

private:
    size_t num_cells_;
    Discrete_Solution_HOM discrete_solution_;
    std::unique_ptr<Governing_Equation> governing_equation_;
};


//#include "Cells_FVM.h"
//#include "Cells_HOM.h"
//#include "Spatial_Discrete_Method.h"
//
//
//template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method>
//class Cells;
//
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Cells<Governing_Equation, FVM, Reconstruction_Method> : public Cells_FVM<Governing_Equation>
//{   
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//
//public:
//    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_FVM<Governing_Equation>(grid) {};
//};
//
//template <typename Governing_Equation, typename Reconstruction_Method>
//class Cells<Governing_Equation, HOM, Reconstruction_Method> : public Cells_HOM<Governing_Equation, Reconstruction_Method>
//{
//private:
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//
//public:
//    Cells(const Grid<space_dimension_>& grid, const Reconstruction_Method& reconstruction_method) : Cells_HOM<Governing_Equation, Reconstruction_Method>(grid, reconstruction_method) {};
//};