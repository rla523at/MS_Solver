#include "../INC/Boundary_Flux_Function_Impl.h"

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

Euclidean_Vector Slip_Wall_BC::calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal)
{
	Euclidean_Vector sol(oc_solution.size());

	const auto p = oc_solution[this->pressure_index_];

	for (ushort i = 0; i < this->space_dimension_; ++i)
	{
		const auto soution_index = 1 + i;
		sol[soution_index] = p * normal[i];
	}

	return sol;
}

void Slip_Wall_BC::calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal)
{
	const auto p = oc_solution[this->pressure_index_];

	for (ushort i = 0; i < this->space_dimension_; ++i)
	{
		bdry_flux_ptr[1 + i] = p * normal[i];
	}
}