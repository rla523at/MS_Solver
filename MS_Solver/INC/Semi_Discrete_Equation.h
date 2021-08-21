#pragma once
#include "Boundaries.h"
#include "Cells.h"
#include "Inner_Faces.h"
#include "Periodic_Boundaries.h"
#include "Pressure_Fix.h"
#include "Numerical_Flux_Function.h"
#include "Post_Solution_Data.h"


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
class Semi_Discrete_Equation
{
private:
    static_require(ms::is_governing_equation<Governing_Equation>, "It should be governing equation");
    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>, "It should be spatial discrete method");
    static_require(ms::is_reconsturction_method<Reconstruction_Method>, "It should be reconstruction method");
    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "It should be numerical flux function");

    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
    static constexpr size_t num_equation_ = Governing_Equation::num_equation();

    using Boundaries_ = Boundaries<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Cells_ = Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
    using Periodic_Boundaries_ = Periodic_Boundaries<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
    using Inner_Faces_ = Inner_Faces<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;

    using Discretized_Solution_ = Cells_::Discretized_Solution_;
    using Residual_ = Cells_::Discretized_Solution_;

private:    
    Reconstruction_Method reconstruction_method_;
    Boundaries_ boundaries_;
    Cells_ cells_;
    Periodic_Boundaries_ periodic_boundaries_;
    Inner_Faces_ inner_faces_;

public:
    Semi_Discrete_Equation(Grid<space_dimension_>&& grid) : reconstruction_method_(std::move(grid)), boundaries_(std::move(grid), reconstruction_method_),
        cells_(grid, reconstruction_method_), periodic_boundaries_(std::move(grid), reconstruction_method_), inner_faces_(std::move(grid), reconstruction_method_) {

        if constexpr (std::is_same_v<Spatial_Discrete_Method, HOM>)
            Post_Solution_Data::initialize_HOM(grid, this->reconstruction_method_);

        if constexpr (ms::can_use_pressure_fix<Governing_Equation, Spatial_Discrete_Method>) {
            SET_TIME_POINT;

            this->boundaries_.initialize_pressure_fix();
            this->cells_.initialize_pressure_fix();
            this->periodic_boundaries_.initialize_pressure_fix();
            this->inner_faces_.initialize_pressure_fix();

            Log::content_ << std::left << std::setw(50) << "@ Pressure Fix precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
            Log::print();
        }

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    };

    template <typename Time_Step_Method>
    double calculate_time_step(const std::vector<Discretized_Solution_>& solutions) const {
        static constexpr double time_step_constant_ = Time_Step_Method::constant();

        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>)
            return this->cells_.calculate_time_step(solutions, time_step_constant_);
        else
            return time_step_constant_;
    }

    //Type1
    std::vector<Residual_> calculate_RHS(const std::vector<Discretized_Solution_>& solutions) {
        const auto num_solution = solutions.size();
        std::vector<Residual_> RHS(num_solution);

        this->boundaries_.calculate_RHS(RHS, solutions);
        this->periodic_boundaries_.calculate_RHS(RHS, solutions);
        this->inner_faces_.calculate_RHS(RHS, solutions);
        this->cells_.calculate_RHS(RHS, solutions);

        return RHS;
    }

    template <typename Initial_Condition>
    std::vector<Discretized_Solution_> calculate_initial_solutions(void) const {
        return this->cells_.calculate_initial_solutions<Initial_Condition>();
    }

    template <typename Initial_Condition>
    void estimate_error(const std::vector<Discretized_Solution_>& computed_solutions, const double time) const {
        this->cells_.estimate_error<Initial_Condition>(computed_solutions, time);
    }

    void reconstruct(std::vector<Discretized_Solution_>& solutions) {
        if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
            this->reconstruction_method_.reconstruct(solutions);

        //if constexpr (ms::can_use_pressure_fix<Governing_Equation, Spatial_Discrete_Method>)
        //    Solution_Scaler::inspect_and_scale(solutions);
    }

};