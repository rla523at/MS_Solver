#pragma once
#include "Euclidean_Vector.h"

class Boundary_Flux_Function
{
public://Query
	virtual Euclidean_Vector calculate(const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) abstract;
	virtual void calculate(double* bdry_flux_ptr, const Euclidean_Vector& oc_solution, const Euclidean_Vector& normal) abstract;
};