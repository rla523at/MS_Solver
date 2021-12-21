#pragma once
#include "Boundaries.h"
#include "Cells.h"
#include "Error.h"
#include "Inner_Faces.h"
#include "Post_Processor.h"
#include "Reconstruction.h"


class Semi_Discrete_Equation
{
public://Command
	virtual Euclidean_Vector_Wrapper discrete_solution_vector_wrapper(void) abstract;
	virtual void reconstruct(void) abstract;

public://Query
	virtual double calculate_time_step(void) const abstract;
	virtual Euclidean_Vector calculate_RHS(void) const abstract;
	virtual Euclidean_Vector discrete_solution_vector(void) const abstract;
	virtual Constant_Euclidean_Vector_Wrapper discrete_solution_constant_vector_wrapper(void) const abstract;
	virtual std::vector<double> calculate_error_norms(const Grid& grid, const double end_time) const abstract;

};

class Semi_Discrete_Equation_DG : public Semi_Discrete_Equation
{
public:
	Semi_Discrete_Equation_DG(const Configuration& configuration, const Grid& grid);
	Semi_Discrete_Equation_DG(std::unique_ptr<Discrete_Solution_DG>&& discrete_solution, std::unique_ptr<Cells_DG>&& cells
		, std::unique_ptr<Boundaries_DG>&& boundaries, std::unique_ptr<Inner_Faces_DG>&& inner_faces)
		: discrete_solution_(std::move(discrete_solution))
		, cells_(std::move(cells))
		, boundaries_(std::move(boundaries))
		, inner_faces_(std::move(inner_faces)) {};

public://Command
	Euclidean_Vector_Wrapper discrete_solution_vector_wrapper(void) override;
	void reconstruct(void) override;

public://Query
	double calculate_time_step(void) const override;
	Euclidean_Vector calculate_RHS(void) const override;
	Euclidean_Vector discrete_solution_vector(void) const override;
	Constant_Euclidean_Vector_Wrapper discrete_solution_constant_vector_wrapper(void) const override;
	std::vector<double> calculate_error_norms(const Grid& grid, const double end_time) const override;

private:
	std::unique_ptr<Discrete_Solution_DG> discrete_solution_;
	std::unique_ptr<Cells_DG> cells_;	
	std::unique_ptr<Boundaries_DG> boundaries_;
	std::unique_ptr<Inner_Faces_DG> inner_faces_;//inter cell faces
	std::unique_ptr<Reconstruction_DG> reconstruction_;
	std::unique_ptr<Error> error_;
};

class Semi_Discrete_Equation_Factory
{
public:
	static std::unique_ptr<Semi_Discrete_Equation> make_unique(const Configuration& configuration, const Grid& grid);
};
































//class Semi_Discrete_Equation
//{
//public:
//	Semi_Discrete_Equation(const Configuration& configuration, const Grid& grid);
//
//public://Command
//	void update_solution(Euclidean_Vector&& updated_soltuion);
//
//public://Query
//	double calculate_time_step(void) const;
//	Euclidean_Vector calculate_RHS(void) const;
//	const Euclidean_Vector& get_solution_vector(void) const;	
//
//private:
//	std::unique_ptr<Cells> cells_;
//	std::vector<std::unique_ptr<Face>> faces_;
//};

























//
//#include "Boundaries.h"
//#include "Cells.h"
//#include "Inner_Faces.h"
//#include "Periodic_Boundaries.h"
//#include "Solution_Scaler.h"
//#include "Numerical_Flux_Function.h"
//#include "Tecplot.h"
//
//
//template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, typename Numerical_Flux_Function, bool scailing_method_flag>
//class Semi_Discrete_Equation
//{
//private:
//    static_require(ms::is_governing_equation<Governing_Equation>, "It should be governing equation");
//    static_require(ms::is_spatial_discrete_method<Spatial_Discrete_Method>, "It should be spatial discrete method");
//    static_require(ms::is_reconsturction_method<Reconstruction_Method>, "It should be reconstruction method");
//    static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "It should be numerical flux function");
//
//    static constexpr size_t space_dimension_ = Governing_Equation::space_dimension();
//    static constexpr size_t num_equation_ = Governing_Equation::num_equation();
//
//    using Boundaries_ = Boundaries<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
//    using Cells_ = Cells<Governing_Equation, Spatial_Discrete_Method, Reconstruction_Method>;
//    using Periodic_Boundaries_ = Periodic_Boundaries<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
//    using Inner_Faces_ = Inner_Faces<Spatial_Discrete_Method, Reconstruction_Method, Numerical_Flux_Function>;
//
//    using Discretized_Solution_ = Cells_::Discretized_Solution_;
//    using Residual_ = Cells_::Discretized_Solution_;
//
//private:    
//    Reconstruction_Method reconstruction_method_;
//    Boundaries_ boundaries_;
//    Cells_ cells_;
//    Periodic_Boundaries_ periodic_boundaries_;
//    Inner_Faces_ inner_faces_;
//
//public:
//    Semi_Discrete_Equation(const Grid<space_dimension_>& grid) : reconstruction_method_(grid), boundaries_(grid, reconstruction_method_),
//        cells_(grid, reconstruction_method_), periodic_boundaries_(grid, reconstruction_method_), inner_faces_(grid, reconstruction_method_) {
//
//        if constexpr (std::is_same_v<Spatial_Discrete_Method, HOM>)
//            Tecplot::initialize_HOM(grid, this->reconstruction_method_); //post
//
//        if constexpr (ms::can_use_scaliling_method<Governing_Equation, Spatial_Discrete_Method>) {
//            SET_TIME_POINT;
//
//            this->boundaries_.initialize_scaling_method();
//            this->cells_.initialize_scaling_method();
//            this->periodic_boundaries_.initialize_scaling_method();
//            this->inner_faces_.initialize_scaling_method();
//
//            Log::content_ << std::left << std::setw(50) << "@ Scailing Method precalculation" << " ----------- " << GET_TIME_DURATION << "s\n\n";
//            Log::print();
//        }
//
//        Log::content_ << "================================================================================\n";
//        Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
//        Log::content_ << "================================================================================\n\n";
//        Log::print();
//    };
//
//    template <typename Time_Step_Method>
//    double calculate_time_step(const std::vector<Discretized_Solution_>& solutions) const {
//        static constexpr double time_step_constant_ = Time_Step_Method::constant();
//
//        if constexpr (std::is_same_v<Time_Step_Method, CFL<time_step_constant_>>)            
//            return this->cells_.calculate_time_step(solutions, time_step_constant_);
//        else
//            return time_step_constant_;
//    }
//        
//    std::vector<Residual_> calculate_RHS(const std::vector<Discretized_Solution_>& solutions) {
//        const auto num_solution = solutions.size();
//        std::vector<Residual_> RHS(num_solution);
//
//        this->boundaries_.calculate_RHS(RHS, solutions);
//        this->periodic_boundaries_.calculate_RHS(RHS, solutions);
//        this->inner_faces_.calculate_RHS(RHS, solutions);
//        this->cells_.calculate_RHS(RHS, solutions);
//
//        return RHS;
//    }
//
//    template <typename Initial_Condition>
//    std::vector<Discretized_Solution_> calculate_initial_solutions(void) const {
//        return this->cells_.calculate_initial_solutions<Initial_Condition>();
//    }
//
//    template <typename Initial_Condition>
//    void estimate_error(const std::vector<Discretized_Solution_>& computed_solutions, const double time) const {
//        if constexpr (ms::is_sine_wave<Initial_Condition> && ms::is_linear_advection<Governing_Equation>)
//            this->cells_.estimate_error(computed_solutions, time);
//        else
//            Log::content_ << "In this setting, error analysis result does not provided.";
//                
//        Log::print();
//    }
//
//    void reconstruct(std::vector<Discretized_Solution_>& solutions) {
//        if constexpr (!ms::is_default_reconstruction<Spatial_Discrete_Method, Reconstruction_Method>)
//            this->reconstruction_method_.reconstruct(solutions);
//
//        if constexpr (scailing_method_flag && ms::can_use_scaliling_method<Governing_Equation, Spatial_Discrete_Method>)
//            Solution_Scaler<space_dimension_>::inspect_and_scale(solutions);
//    }
//
//};