#include "../INC/Semi_Discrete_Equation.h"

Semi_Discrete_Equation_DG::Semi_Discrete_Equation_DG(const Configuration& configuration, const Grid& grid) :
	discrete_solution_(configuration, grid),
	cells_(configuration, grid, this->discrete_solution_) {};

void Semi_Discrete_Equation_DG::update_solution(Euclidean_Vector&& updated_soltuion_v)
{
	this->discrete_solution_.update_solution(std::move(updated_soltuion_v));
}

double Semi_Discrete_Equation_DG::calculate_time_step(void) const
{
	return this->cells_.calculate_time_step(this->discrete_solution_);
}

Euclidean_Vector Semi_Discrete_Equation_DG::calculate_RHS(void) const
{
	std::vector<double> RHS(this->discrete_solution_.num_values());

	this->cells_.calculate_RHS(RHS.data(),this->discrete_solution_);
}

const Euclidean_Vector& Semi_Discrete_Equation_DG::get_solution_vector(void) const
{
	return this->discrete_solution_.get_solution_vector();
}

std::unique_ptr<Semi_Discrete_Equation> Semi_Discrete_Equation_Factory::make(const Configuration& configuration, const Grid& grid)
{
	const auto name = configuration.get("spatial_discrete_scheme");

	if (ms::contains_icase(name, "DG"))
	{
		return std::make_unique<Semi_Discrete_Equation_DG>(configuration, grid);
	}
	else
	{
		EXCEPTION("spatial discrete scheme in configuration is not supported");
	}
}