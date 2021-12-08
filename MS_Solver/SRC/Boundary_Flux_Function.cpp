#include "../INC/Boundary_Flux_Function.h"

Euclidean_Vector Initial_Constant_BC::calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal)
{
	this->calculate_neighbor_solution(oc_solution);
	return this->numerical_flux_function_->calculate(oc_solution, this->neighbor_solution_, normal);
}

void Initial_Constant_BC::calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal)
{
	this->calculate_neighbor_solution(oc_solution);
	this->numerical_flux_function_->calculate(bdry_flux_ptr, oc_solution, this->neighbor_solution_, normal);
}

void Initial_Constant_BC::calculate_neighbor_solution(const Euclidean_Vector& oc_solution)
{
	if (!this->is_initialized_)
	{
		this->neighbor_solution_ = oc_solution;
		this->is_initialized_ = true;
	}
}

std::unique_ptr<Boundary_Flux_Function> Boundary_Flux_Function_Factory::make_unique(const ElementType boundary_type, const std::shared_ptr<Numerical_Flux_Function>& numerical_flux_function)
{
	switch (boundary_type)
	{
	case ElementType::initial_constant_BC:
	{
		return std::make_unique<Initial_Constant_BC>(numerical_flux_function);
	}
	default:
	{
		EXCEPTION("not supported boundary type");
		return nullptr;
	}
	}
};