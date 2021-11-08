#include "../INC/Semi_Discrete_Equation.h"

void Semi_Discrete_Equation::update_solution(Euclidean_Vector&& updated_soltuion_v)
{
	this->cells_->update_solution(std::move(updated_soltuion_v));
}

double Semi_Discrete_Equation::calculate_time_step(void) const
{
	return this->cells_->calculate_time_step();
}

Euclidean_Vector Semi_Discrete_Equation::calculate_RHS(void) const
{
	const auto num_solution_values = this->cells_->num_solution_values();
	std::vector<double> RHS(num_solution_values);

	this->cells_->calculate_RHS(RHS.data());
}

const Euclidean_Vector& Semi_Discrete_Equation::get_solution_vector(void) const
{
	return this->cells_->get_solution_vector();
}
