#pragma once
#include "Euclidean_Vector.h"

#include <numbers>

class Initial_Condition
{
public:
	virtual Euclidean_Vector calculate_solution(const Euclidean_Vector& space_vector) const abstract;
	virtual double target_end_time(void) const { return -1.0; };
};