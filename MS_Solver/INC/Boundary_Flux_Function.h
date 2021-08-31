#pragma once
#include "Numerical_Flux_Function.h"
#include "Element.h"

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
class Supersonic_Inlet : public Boundary_Flux_Function<Numerical_Flux_Function>
{
private:
	using This_			= Supersonic_Inlet<Numerical_Flux_Function>;
	using Space_Vector_ = This_::Space_Vector_;
	using Solution_		= This_::Solution_;

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
class Supersonic_Outlet : public Boundary_Flux_Function<Numerical_Flux_Function>
{
private:
	using This_			= Supersonic_Outlet<Numerical_Flux_Function>;
	using Space_Vector_ = This_::Space_Vector_;
	using Solution_		= This_::Solution_;

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


template <typename Numerical_Flux_Function>
std::unique_ptr<Boundary_Flux_Function<Numerical_Flux_Function>> Boundary_Flux_Function_Factory<Numerical_Flux_Function>::make(const ElementType boundary_type) {
	switch (boundary_type)
	{
	case ElementType::supersonic_inlet:
		return std::make_unique<Supersonic_Inlet<Numerical_Flux_Function>>();
	case ElementType::supersonic_outlet:
		return std::make_unique<Supersonic_Outlet<Numerical_Flux_Function>>();
	//case ElementType::slip_wall_2D:
	//	return std::make_unique<Slip_Wall_2D>();
	default:
		throw std::runtime_error("wrong element type");
		break;
	}	
}