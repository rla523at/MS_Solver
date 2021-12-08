#pragma once
#include "Euclidean_Vector.h"
#include "Element.h"
#include "Numerical_Flux_Function.h"

class Boundary_Flux_Function
{
public://Query
	virtual Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) abstract;
	virtual void calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) abstract;

};

class Initial_Constant_BC : public Boundary_Flux_Function
{
public:
	Initial_Constant_BC(const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
		: numerical_flux_function_(numerical_flux_function) {};

public://Query
	Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) override;
	void calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) override;

private:
	void calculate_neighbor_solution(const Euclidean_Vector& oc_solution);


private:
	bool is_initialized_ = false;
	Euclidean_Vector neighbor_solution_;
	std::shared_ptr<Numerical_Flux_Function> numerical_flux_function_;
};


class Boundary_Flux_Function_Factory//static class
{
public:
	static std::unique_ptr<Boundary_Flux_Function> make_unique(const ElementType boundary_type, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function);

private:
	Boundary_Flux_Function_Factory(void) = delete;
};






//template <ushort num_equation>
//class Supersonic_Inlet1_Neighbor_Solution_Calculator
//{
//private:
//	Supersonic_Inlet1_Neighbor_Solution_Calculator(void) = delete;
//
//private:
//	inline static Euclidean_Vector<num_equation> inflow_;
//
//public:
//	static void initialize(const Euclidean_Vector<num_equation>& inflow) { inflow_ = inflow; };
//	static Euclidean_Vector<num_equation> calculate(void) { return inflow_; };
//};
//
//
//template <ushort num_equation>
//class Supersonic_Inlet2_Neighbor_Solution_Calculator
//{
//private:
//	Supersonic_Inlet2_Neighbor_Solution_Calculator(void) = delete;
//
//private:
//	inline static Euclidean_Vector<num_equation> inflow_;
//
//public:
//	static void initialize(const Euclidean_Vector<num_equation>& inflow) { inflow_ = inflow; };
//	static Euclidean_Vector<num_equation> calculate(void) { return inflow_; };
//};

//template <ushort num_equation>
//class Slip_Wall_Neighbor_Solution_Calculator
//{
//private:
//	Slip_Wall_Neighbor_Solution_Calculator(void) = delete;
//
//public:
//	static Euclidean_Vector<num_equation> calculate(const Euclidean_Vector<num_equation>& solution) {
//		std::array<double, num_equation> values = -1 * solution;
//		values.front() *= -1;
//		values.back() *= -1;
//		return values;
//	};
//};
//
//
//template <ushort num_equation>
//class Initial_Constant_BC_Neighbor_Solution_Calculator
//{
//private:
//	bool is_first_ = true;
//	Euclidean_Vector<num_equation> initial_constant_;
//
//public:
//	Euclidean_Vector<num_equation> calculate(const Euclidean_Vector<num_equation>& solution) {
//		if (this->is_first_) {
//			this->initial_constant_ = solution;
//			this->is_first_ = false;
//		}
//
//		return this->initial_constant_;
//	};
//};



//template <typename Numerical_Flux_Function>
//class Supersonic_Inlet1 : public Boundary_Flux_Function<Numerical_Flux_Function>
//{
//private:
//	using This_			= Supersonic_Inlet1<Numerical_Flux_Function>;
//	using Space_Vector_ = This_::Space_Vector_;
//	using Solution_		= This_::Solution_;
//
//public:
//	This_::Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) const override {
//		return Numerical_Flux_Function::calculate(solution, Supersonic_Inlet1_Neighbor_Solution_Calculator<This_::num_equation_>::calculate(), normal);
//	}
//};


//template <typename Numerical_Flux_Function>
//class Supersonic_Inlet2 : public Boundary_Flux_Function<Numerical_Flux_Function>
//{
//private:
//	using This_			= Supersonic_Inlet2<Numerical_Flux_Function>;
//	using Space_Vector_ = This_::Space_Vector_;
//	using Solution_		= This_::Solution_;
//
//public:
//	This_::Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) const override {
//		return Numerical_Flux_Function::calculate(solution, Supersonic_Inlet2_Neighbor_Solution_Calculator<This_::num_equation_>::calculate(), normal);
//	}
//};


//template <typename Numerical_Flux_Function>
//class Supersonic_Outlet : public Boundary_Flux_Function<Numerical_Flux_Function>
//{
//private:
//	using This_			= Supersonic_Outlet<Numerical_Flux_Function>;
//	using Space_Vector_ = This_::Space_Vector_;
//	using Solution_		= This_::Solution_;
//
//public:
//	This_::Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) const override {
//		return Numerical_Flux_Function::calculate(solution, solution, normal);
//	}
//};


//template <typename Numerical_Flux_Function>
//class Slip_Wall : public Boundary_Flux_Function<Numerical_Flux_Function>
//{
//private:
//	using This_			= Slip_Wall<Numerical_Flux_Function>;
//	using Space_Vector_ = This_::Space_Vector_;
//	using Solution_		= This_::Solution_;
//
//public:
//	This_::Boundary_Flux_ calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) const override {
//		static_require(This_::space_dimension_ <= 3, "size can not exceed 3");
//
//		if constexpr (ms::is_Euler<Numerical_Flux_Function::Governing_Equation_>) {
//			const auto pvariable = Euler<This_::space_dimension_>::conservative_to_primitive(oc_cvariable);
//
//			if constexpr (This_::space_dimension_ == 2) {
//				const auto p = pvariable[2];
//				return { 0.0, p * normal[0], p * normal[1], 0.0 };
//			}
//			else {
//				const auto p = pvariable[3];
//				return { 0.0, p * normal[0], p * normal[1], p * normal[2], 0.0 };
//			}
//		}
//		else {
//			throw std::runtime_error("Governing equation should be Euler");			
//			return {};
//		}
//	}
//};




