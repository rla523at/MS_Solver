#pragma once
#include "Boundaries.h"
#include "Cells.h"
#include "Inner_Faces.h"
#include "Periodic_Boundaries.h"
#include "Pressure_Fix.h"
#include "Numerical_Flux_Function.h"


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

private:
    Semi_Discrete_Equation(void) = delete;

public:
    static void initialize(Grid<space_dimension_>&& grid) 
        : reconstruction_method_(grid), boundaries_(std::move(grid), reconstruction_method_), cells_(grid, reconstruction_method_){

        Cells_::initialize(std::move(grid));
        Boundaries_::initialize(std::move(grid));
        Periodic_Boundaries_::initialize(std::move(grid));
        Inner_Faces_::initialize(std::move(grid));

        if constexpr (ms::can_use_pressure_fix<Governing_Equation, Spatial_Discrete_Method>)
            Pressure_Fix::initialize<Reconstruction_Method>(grid);

        Log::content_ << "================================================================================\n";
        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        Log::content_ << "================================================================================\n\n";
        Log::print();
    };

    template <typename Time_Step_Method>
    static double calculate_time_step(const std::vector<Discretized_Solution_>& solutions) {
        static constexpr double time_step_constant_ = Time_Step_Method::constant();
        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>)
            return Cells_::calculate_time_step(solutions, time_step_constant_);
        else
            return time_step_constant_;
    }

    //Type1
    static std::vector<Residual_> calculate_RHS(std::vector<Discretized_Solution_>& solutions) {
        const auto num_solution = solutions.size();
        std::vector<Residual_> RHS(num_solution);

        if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
            reconstruction_method_.reconstruct(solutions);

        Boundaries_::calculate_RHS(RHS, solutions);
        Periodic_Boundaries_::calculate_RHS(RHS, solutions);
        Inner_Faces_::calculate_RHS(RHS, solutions);
        Cells_::calculate_RHS(RHS, solutions);

        return RHS;
    }

    template <typename Initial_Condition>
    static std::vector<Discretized_Solution_> calculate_initial_solutions(void) {
        return Cells_::template calculate_initial_solutions<Initial_Condition>();
    }

    template <typename Initial_Condition>
    static void estimate_error(const std::vector<Discretized_Solution_>& computed_solutions, const double time) {
        Cells_::template estimate_error<Initial_Condition>(computed_solutions, time);
    }
};































//template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function>
//class Semi_Discrete_Equation
//{
//private:
//    static_require(ms::is_governing_equation<Governing_Equation>,               "It should be governing equation");
//    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>,     "It should be spatial discrete method");
//    static_require(ms::is_reconsturction_method<Reconstruction_Method>,         "It should be reconstruction method");
//    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>,    "It should be numerical flux function");
//
//    static constexpr size_t space_dimension_    = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_       = Governing_Equation::num_equation();    
//
//    using Boundaries_           = Boundaries<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
//    using Cells_                = Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
//    using Periodic_Boundaries_  = Periodic_Boundaries<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
//    using Inner_Faces_          = Inner_Faces<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
//
//    using Discretized_Solution_ = Cells_::Discretized_Solution_;
//    using Residual_             = Cells_::Discretized_Solution_;
//
//private:
//    Reconstruction_Method reconstruction_method_;
//
//private:
//    Semi_Discrete_Equation(void) = delete;
//    
//public:
//    static void initialize(Grid<space_dimension_>&& grid) :reconstruction_method_(grid) {
//        /*Reconstruction_Method::initialize(std::move(grid));*/
//        Cells_::initialize(std::move(grid));
//        Boundaries_::initialize(std::move(grid));
//        Periodic_Boundaries_::initialize(std::move(grid));
//        Inner_Faces_::initialize(std::move(grid));
//
//        if constexpr (ms::can_use_pressure_fix<Governing_Equation, Spatial_Discrete_Method>)
//            Pressure_Fix::initialize<Reconstruction_Method>(grid);
//
//        Log::content_ << "================================================================================\n";
//        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
//        Log::content_ << "================================================================================\n\n";
//        Log::print();
//    };
//
//    template <typename Time_Step_Method>
//    static double calculate_time_step(const std::vector<Discretized_Solution_>& solutions) {
//        static constexpr double time_step_constant_ = Time_Step_Method::constant();
//        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>) 
//            return Cells_::calculate_time_step(solutions, time_step_constant_);
//        else
//            return time_step_constant_;
//    }
//
//    //Type1
//    static std::vector<Residual_> calculate_RHS(std::vector<Discretized_Solution_>& solutions) {
//        const auto num_solution = solutions.size();
//        std::vector<Residual_> RHS(num_solution);
//
//        if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
//            Reconstruction_Method::reconstruct(solutions);
//
//        //Pressure_Fix::reconstruct(solutions);
//
//        Boundaries_::calculate_RHS(RHS, solutions);
//        Periodic_Boundaries_::calculate_RHS(RHS, solutions);
//        Inner_Faces_::calculate_RHS(RHS, solutions);
//        Cells_::calculate_RHS(RHS, solutions);
//
//        return RHS;
//    }
//
//    ////Type2
//    //static std::vector<Residual_> calculate_RHS(const std::vector<Discretized_Solution_>& solutions) {
//    //    static const auto num_solution = solutions.size();
//    //    std::vector<Residual_> RHS(num_solution);
//
//    //    //if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
//    //    auto reconstructed_solution = solutions;
//    //    Reconstruction_Method::reconstruct(reconstructed_solution);
//
//    //    Boundaries_::calculate_RHS(RHS, reconstructed_solution);
//    //    Periodic_Boundaries_::calculate_RHS(RHS, reconstructed_solution);
//    //    Inner_Faces_::calculate_RHS(RHS, reconstructed_solution);
//    //    Cells_::calculate_RHS(RHS, reconstructed_solution);
//
//    //    return RHS;
//    //}
//
//    //static std::vector<Residual_> calculate_RHS(const std::vector<Discretized_Solution_>& solutions) {
//    //    static const auto num_solution = solutions.size();
//    //    std::vector<Residual_> RHS(num_solution);
//
//    //    if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
//    //        Reconstruction_Method::reconstruct(solutions);
//
//    //    Boundaries_::calculate_RHS(RHS, solutions);
//    //    Periodic_Boundaries_::calculate_RHS(RHS, solutions);
//    //    Inner_Faces_::calculate_RHS(RHS, solutions);
//    //    Cells_::calculate_RHS(RHS, solutions);
//
//    //    return RHS;
//    //}
//
//    template <typename Initial_Condition>
//    static std::vector<Discretized_Solution_> calculate_initial_solutions(void) {
//        return Cells_::template calculate_initial_solutions<Initial_Condition>();
//    }
//
//    template <typename Initial_Condition>
//    static void estimate_error(const std::vector<Discretized_Solution_>& computed_solutions, const double time) {
//        Cells_::template estimate_error<Initial_Condition>(computed_solutions, time);
//    }
//};
//
