#pragma once
#include "Boundaries.h"
#include "Cells.h"
#include "Inner_Faces.h"
#include "Periodic_Boundaries.h"
#include "Numerical_Flux_Function.h"

template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Semi_Discrete_Equation
{
private:
    static_require(ms::is_governing_equation<Governing_Equation>,               "It should be governing equation");
    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>,     "It should be spatial discrete method");
    static_require(ms::is_reconsturction_method<Reconstruction_Method>,         "It should be reconstruction method");
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>,    "It should be numerical flux function");

    static constexpr size_t space_dimension_    = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_       = Governing_Equation::num_equation();

    using Boundaries_           = Boundaries<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Cells_                = Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Periodic_Boundaries_  = Periodic_Boundaries<Spatial_Discrete_Method, Reconstruction_Method, space_dimension_>;
    using Inner_Faces_          = Inner_Faces<Spatial_Discrete_Method, Reconstruction_Method, space_dimension_>;

    using Solution_             = typename Governing_Equation::Solution_;
    using Residual_             = Euclidean_Vector<num_equation_>;

private:
    Semi_Discrete_Equation(void) = delete;
    
public:
    static void initialize(Grid<space_dimension_>&& grid) {
        Reconstruction_Method::initialize(std::move(grid));
        Cells_::initialize(std::move(grid));
        Boundaries_::initialize(std::move(grid));
        Periodic_Boundaries_::initialize(std::move(grid));
        Inner_Faces_::initialize(std::move(grid));


        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    };

    template <typename Time_Step_Method>
    static double calculate_time_step(const std::vector<Solution_>& solutions) {
        static constexpr double time_step_constant_ = Time_Step_Method::constant();
        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>) {
            const auto projected_maximum_lambdas = Governing_Equation::coordinate_projected_maximum_lambdas(solutions);
            return Cells_::calculate_time_step(projected_maximum_lambdas, time_step_constant_);
        }
        else
            return time_step_constant_;
    }

    static std::vector<Residual_> calculate_RHS(const std::vector<Solution_>& solutions) {
        static const auto num_solution = solutions.size();
        std::vector<Residual_> RHS(num_solution);

        if constexpr (!ms::is_default_reconstruction<Reconstruction_Method>)
            Reconstruction_Method::reconstruct(solutions);

        Boundaries_::calculate_RHS(RHS, solutions);
        Periodic_Boundaries_::template calculate_RHS<Numerical_Flux_Function>(RHS, solutions);
        Inner_Faces_::template calculate_RHS<Numerical_Flux_Function>(RHS, solutions);
        Cells_::calculate_RHS(RHS, solutions);

        return RHS;
    }

    template <typename Initial_Condition>
    static std::vector<Solution_> calculate_initial_solutions(void) {
        return Cells_::template calculate_initial_solutions<Initial_Condition>();
    }

    template <typename Initial_Condition>
    static void estimate_error(const std::vector<Solution_>& computed_solution, const double time) {
        Cells_::template estimate_error<Initial_Condition, Governing_Equation>(computed_solution, time);
    }

};