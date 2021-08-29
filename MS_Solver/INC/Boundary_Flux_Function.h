#pragma once
#include "Numerical_Flux_Function.h"
#include "Element.h"
#include "Setting.h"

template <typename Numerical_Flux_Function>
class Boundary_Flux_Function
{
protected:
	static_require(ms::is_numeirical_flux_function<Numerical_Flux_Function>, "it shuld be numerical flux function");

	static constexpr ushort space_dimension_	= Numerical_Flux_Function::space_dimension();
	static constexpr ushort num_equation_		= Numerical_Flux_Function::num_equation();

	using Space_Vector_		= Euclidean_Vector<space_dimension_>;
	using Solution_			= Euclidean_Vector<num_equation_>;
	using Boundary_Flux_	= Euclidean_Vector<num_equation_>;

public:
	virtual Boundary_Flux_ calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) const abstract;
};


template <typename Numerical_Flux_Function>
class Supersonic_Inlet_2D : public Boundary_Flux_Function<Numerical_Flux_Function>
{
private:
	using This_ = Supersonic_Inlet_2D<Numerical_Flux_Function>;
	using Space_Vector_ = This_::Space_Vector_;
	using Solution_ = This_::Solution_;

private:
	inline static This_::Solution_ inflow_;

public:
	static void initialize(const Solution_& inflow) { This_::inflow_ = inflow; };

public:
	This_::Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) const override {
		return Numerical_Flux_Function::calculate(solution, This_::inflow_, normal);
	}
};


template <typename Numerical_Flux_Function>
class Supersonic_Outlet_2D : public Boundary_Flux_Function<Numerical_Flux_Function>
{
private:
	using This_ = Supersonic_Outlet_2D<Numerical_Flux_Function>;
	using Space_Vector_ = This_::Space_Vector_;
	using Solution_ = This_::Solution_;

public:
	This_::Boundary_Flux_ calculate(const Solution_& solution, const Space_Vector_& normal) const override {
		return Numerical_Flux_Function::calculate(solution, solution, normal);
	}

};


//class Slip_Wall_2D : public Boundary_Flux_Function<LLF<Euler_2D>>
//{
//	Boundary_Flux_ calculate(const Solution_& oc_cvariable, const Space_Vector_& normal) const override;
//};



template <typename Numerical_Flux_Function>
class Boundary_Flux_Function_Factory
{
public:
	static std::unique_ptr<Boundary_Flux_Function<Numerical_Flux_Function>> make(const ElementType boundary_type);
};



//template definition part
//template <typename Governing_Equation> 
//Supersonic_Outlet_2D<Governing_Equation>::Boundary_Flux_ Supersonic_Outlet_2D<Governing_Equation>::calculate(const Solution_& solution, const Space_Vector_& normal) const {
//	////reflective boundaries
//	//std::array<double, num_equation_> boundary_values;
//
//	//boundary_values[0] = solution[0];
//	//boundary_values[1] = -solution[1];
//	//boundary_values[2] = -solution[2];
//	//boundary_values[3] = solution[3];
//
//	//return LLF<Euler_2D>::calculate(solution, boundary_values, normal);
//
//
//	//if (is_first) { //debug
//	//	initial_solution = solution;
//	//	is_first = false;
//	//}
//	//
//	//const auto oc_physical_flux = Governing_Equation::physical_flux(initial_solution);
//	//return oc_physical_flux * normal;//debug
//	
//	const auto oc_physical_flux = Governing_Equation::physical_flux(solution);
//	return oc_physical_flux * normal;
//}


template <typename Numerical_Flux_Function>
std::unique_ptr<Boundary_Flux_Function<Numerical_Flux_Function>> Boundary_Flux_Function_Factory<Numerical_Flux_Function>::make(const ElementType boundary_type) {
	switch (boundary_type)
	{
	case ElementType::supersonic_inlet_2D:
		return std::make_unique<Supersonic_Inlet_2D<Numerical_Flux_Function>>();
	case ElementType::supersonic_outlet_2D:
		return std::make_unique<Supersonic_Outlet_2D<Numerical_Flux_Function>>();
	//case ElementType::slip_wall_2D:
	//	return std::make_unique<Slip_Wall_2D>();
	default:
		throw std::runtime_error("wrong element type");
		break;
	}
	
}