#pragma once
#include "Boundary_Flux_Function.h"


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


class Slip_Wall_BC : public Boundary_Flux_Function
{
public:
	Slip_Wall_BC(const std::shared_ptr<Governing_Equation>& governing_equation, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
		: pressure_index_(governing_equation->pressure_index())
		, space_dimension_(governing_equation->space_dimension())
		, numerical_flux_function_(numerical_flux_function) {};

public://Query
	Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) override
	{
		Euclidean_Vector sol(oc_solution.size());

		const auto p = oc_solution[this->pressure_index_];

		for (ushort i = 0; i < this->space_dimension_; ++i)
		{
			sol[1 + i] = p * normal[i];
		}

		return sol;
	}

	void calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) override
	{
		const auto p = oc_solution[this->pressure_index_];

		for (ushort i = 0; i < this->space_dimension_; ++i)
		{
			bdry_flux_ptr[1 + i] = p * normal[i];
		}
	}

private:
	ushort pressure_index_ = 0;
	ushort space_dimension_ = 0;
	std::shared_ptr<Numerical_Flux_Function> numerical_flux_function_;
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

