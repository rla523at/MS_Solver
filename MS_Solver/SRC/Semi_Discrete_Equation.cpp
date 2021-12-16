#include "../INC/Semi_Discrete_Equation.h"

Semi_Discrete_Equation_DG::Semi_Discrete_Equation_DG(const Configuration& configuration, const Grid& grid)
{
	const auto governing_equation = Governing_Equation_Factory::make_shared(configuration);

	const auto initial_condition = Initial_Condition_Factory::make_unique(configuration);
	const auto solution_degree = configuration.solution_degree();
	this->discrete_solution_ = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

	auto time_step_calculator = Time_Step_Calculator_Factory::make_unique(configuration, grid);
	this->cells_ = std::make_unique<Cells_DG>(governing_equation, std::move(time_step_calculator), grid, *this->discrete_solution_);

	const auto numerical_flux_function = Numerical_Flux_Function_Factory::make_shared(configuration, governing_equation);
	this->boundaries_ = std::make_unique<Boundaries_DG>(grid, *this->discrete_solution_, numerical_flux_function);
	this->inner_faces_ = std::make_unique<Inner_Faces_DG>(numerical_flux_function, grid, *this->discrete_solution_);

	this->reconstruction_ = Reconstruction_DG_Factory::make_unique(configuration, grid, *this->discrete_solution_);
	this->error_ = Error_Factory::make_unqiue(configuration);

	Post_Processor::initialize(configuration, grid, *this->discrete_solution_);
}

Euclidean_Vector_Wrapper Semi_Discrete_Equation_DG::discrete_solution_vector_wrapper(void)
{
	return this->discrete_solution_->discrete_solution_vector_wrapper();
}

void Semi_Discrete_Equation_DG::reconstruct(void)
{
	return this->reconstruction_->reconstruct(*this->discrete_solution_);
}

void Semi_Discrete_Equation_DG::update_solution(const Euclidean_Vector& updated_soltuion_v)
{
	this->discrete_solution_->update_solution(updated_soltuion_v);
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
	Residual RHS(this->discrete_solution_->num_total_values(), this->discrete_solution_->get_coefficient_start_indexes());
	this->cells_->calculate_RHS(RHS, *this->discrete_solution_);
	this->boundaries_->calculate_RHS(RHS, *this->discrete_solution_);
	this->inner_faces_->calculate_RHS(RHS, *this->discrete_solution_);

	return RHS.move_values();	//���׸� ����ų���� �����Ͱ�����
}

Euclidean_Vector Semi_Discrete_Equation_DG::discrete_solution_vector(void) const
{
	return this->discrete_solution_->discrete_solution_vector();
}

Euclidean_Vector_Constant_Wrapper Semi_Discrete_Equation_DG::solution_vector_constant_wrapper(void) const
{
	return this->discrete_solution_->solution_vector_constant_wrapper();
}

std::vector<double> Semi_Discrete_Equation_DG::calculate_error_norms(const Grid& grid, const double end_time) const
{
	return this->error_->calculate_error_norms(grid, *this->discrete_solution_, end_time);
}

std::unique_ptr<Semi_Discrete_Equation> Semi_Discrete_Equation_Factory::make_unique(const Configuration& configuration, const Grid& grid)
{
	const auto& name = configuration.get_spatial_discrete_scheme();

	if (ms::contains_icase(name, "DG"))
	{
		return std::make_unique<Semi_Discrete_Equation_DG>(configuration, grid);

		//const auto governing_equation = Governing_Equation_Factory::make_shared(configuration);

		//const auto initial_condition = Initial_Condition_Factory::make_unique(configuration);
		//const auto solution_degree = configuration.get<ushort>("solution_degree");
		//auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		//auto time_step_calculator = Time_Step_Calculator_Factory::make_unique(configuration, grid);
		//auto cells_DG = std::make_unique<Cells_DG>(governing_equation, std::move(time_step_calculator), grid, *discrete_solution_DG);

		//const auto numerical_flux_function = Numerical_Flux_Function_Factory::make_shared(configuration, governing_equation);
		//auto boundaries_DG = std::make_unique<Boundaries_DG>(grid, *discrete_solution_DG, numerical_flux_function);
		//auto inner_faces_DG = std::make_unique<Inner_Faces_DG>(numerical_flux_function, grid, *discrete_solution_DG);

		//return std::make_unique<Semi_Discrete_Equation_DG>(std::move(discrete_solution_DG), std::move(cells_DG), std::move(boundaries_DG), std::move(inner_faces_DG));
	}
	else
	{
		EXCEPTION("spatial discrete scheme in configuration is not supported");
		return nullptr;
	}
}