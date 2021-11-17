#include "../INC/Semi_Discrete_Equation.h"

Semi_Discrete_Equation_DG::Semi_Discrete_Equation_DG(const Configuration& configuration, const Grid& grid)
{
	const auto governing_equation = Governing_Equation_Factory::make_shared(configuration);

	const auto initial_condition = Initial_Condition_Factory::make_unique(configuration);
	const auto solution_degree = configuration.get<ushort>("solution_degree");
	this->discrete_solution_ = std::make_unique<Discrete_Solution_DG>(grid, governing_equation, *initial_condition, solution_degree);
	
	auto time_step_calculator = Time_Step_Calculator_Factory::make_unique(configuration, grid);
	this->cells_ = std::make_unique<Cells_DG>(governing_equation, std::move(time_step_calculator), grid, *this->discrete_solution_);

	const auto numerical_flux_function = Numerical_Flux_Function_Factory::make_shared(configuration, governing_equation);
	this->boundaries_ = std::make_unique<Boundaries_DG>(grid, discrete_solution_, numerical_flux_function);
	this->inner_faces_ = std::make_unique<Inner_Faces_DG>(grid, discrete_solution_, numerical_flux_function);
}

void Semi_Discrete_Equation_DG::update_solution(Euclidean_Vector&& updated_soltuion_v)
{
	this->discrete_solution_->update_solution(std::move(updated_soltuion_v));
}

double Semi_Discrete_Equation_DG::calculate_time_step(void) const
{
	return this->cells_->calculate_time_step(*this->discrete_solution_);
}

Euclidean_Vector Semi_Discrete_Equation_DG::calculate_RHS(void) const
{
	Residual RHS(this->discrete_solution_->num_values(), this->discrete_solution_->get_coefficient_start_indexes());
	this->cells_->calculate_RHS(RHS, *this->discrete_solution_);
	this->boundaries_->calculate_RHS(RHS, *this->discrete_solution_);
	this->inner_faces_->calculate_RHS(RHS, *this->discrete_solution_);

	return RHS.move_values();	//버그를 일으킬수도 있을것같은데
}

const Euclidean_Vector& Semi_Discrete_Equation_DG::get_solution_vector(void) const
{
	return this->discrete_solution_->get_solution_vector();
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