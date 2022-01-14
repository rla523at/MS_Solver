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
	this->inner_faces_ = std::make_unique<Inner_Faces_DG>(grid, *this->discrete_solution_, numerical_flux_function);

	auto boundary_flux_functions = Boundary_Flux_Function_Factory::make_bdry_flux_functions(grid, governing_equation, numerical_flux_function);
	this->boundaries_ = std::make_unique<Boundaries_DG>(grid, *this->discrete_solution_, std::move(boundary_flux_functions));

	this->RHS_ = std::make_unique<Residual>(this->discrete_solution_->num_total_values(), this->discrete_solution_->get_coefficient_start_indexes());
		
	this->reconstruction_ = Reconstruction_DG_Factory::make_unique(configuration, grid, *this->discrete_solution_, governing_equation);
	this->error_ = Error_Factory::make_unqiue(configuration);

	Post_Processor::initialize(configuration, grid, *this->discrete_solution_);
}

void Semi_Discrete_Equation_DG::calculate_RHS(void)
{
	this->RHS_->initialize();
	this->cells_->calculate_RHS(*this->RHS_, *this->discrete_solution_);
	this->boundaries_->calculate_RHS(*this->RHS_, *this->discrete_solution_);
	this->inner_faces_->calculate_RHS(*this->RHS_, *this->discrete_solution_);
}

Euclidean_Vector_Wrapper Semi_Discrete_Equation_DG::discrete_solution_vector_wrapper(void)
{
	return this->discrete_solution_->discrete_solution_vector_wrapper();
}

void Semi_Discrete_Equation_DG::reconstruct(void)
{
	return this->reconstruction_->reconstruct(*this->discrete_solution_);
}

double Semi_Discrete_Equation_DG::calculate_time_step(void) const
{
	return this->cells_->calculate_time_step(*this->discrete_solution_);
}


Euclidean_Vector Semi_Discrete_Equation_DG::discrete_solution_vector(void) const
{
	return this->discrete_solution_->discrete_solution_vector();
}

Constant_Euclidean_Vector_Wrapper Semi_Discrete_Equation_DG::discrete_solution_constant_vector_wrapper(void) const
{
	return this->discrete_solution_->discrete_solution_constant_vector_wrapper();
}

Constant_Euclidean_Vector_Wrapper Semi_Discrete_Equation_DG::RHS_constant_vector_wrapper(void) const
{
	return this->RHS_->residual_constant_vector_wrapper();
}

Euclidean_Vector Semi_Discrete_Equation_DG::RHS_vector(void) const
{
	return this->RHS_->residual_vector();
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