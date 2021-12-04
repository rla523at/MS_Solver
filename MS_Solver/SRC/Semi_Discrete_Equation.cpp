#include "../INC/Semi_Discrete_Equation.h"

Semi_Discrete_Equation_DG::Semi_Discrete_Equation_DG(const Configuration& configuration, const Grid& grid)
{
	const auto governing_equation = Governing_Equation_Factory::make_shared(configuration);

	const auto initial_condition = Initial_Condition_Factory::make_unique(configuration);
	const auto solution_degree = configuration.get<ushort>("solution_degree");
	this->discrete_solution_ = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

	auto time_step_calculator = Time_Step_Calculator_Factory::make_unique(configuration, grid);
	this->cells_ = std::make_unique<Cells_DG>(governing_equation, std::move(time_step_calculator), grid, *this->discrete_solution_);

	const auto numerical_flux_function = Numerical_Flux_Function_Factory::make_shared(configuration, governing_equation);
	this->boundaries_ = std::make_unique<Boundaries_DG>(grid, *this->discrete_solution_, numerical_flux_function);
	this->inner_faces_ = std::make_unique<Inner_Faces_DG>(numerical_flux_function, grid, *this->discrete_solution_);

	Post_Processor::initialize(configuration, grid, *this->discrete_solution_);
}

Euclidean_Vector_Wrapper Semi_Discrete_Equation_DG::solution_vector_wrapper(void)
{
	return this->discrete_solution_->solution_vector_wrapper();
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
	Residual RHS(this->discrete_solution_->num_values(), this->discrete_solution_->get_coefficient_start_indexes());
	this->cells_->calculate_RHS(RHS, *this->discrete_solution_);
	this->boundaries_->calculate_RHS(RHS, *this->discrete_solution_);
	this->inner_faces_->calculate_RHS(RHS, *this->discrete_solution_);

	return RHS.move_values();	//���׸� ����ų���� �����Ͱ�����
}

Euclidean_Vector Semi_Discrete_Equation_DG::solution_vector(void) const
{
	return this->discrete_solution_->solution_vector();
}

Euclidean_Vector_Constant_Wrapper Semi_Discrete_Equation_DG::solution_vector_constant_wrapper(void) const
{
	return this->discrete_solution_->solution_vector_constant_wrapper();
}

std::vector<double> Semi_Discrete_Equation_DG::calculate_error_values(const Exact_Solution& exact_solution, const Grid& grid, const double end_time) const
{
	const auto num_cells = grid.num_cells();

	auto arithmetic_mean_L1_error = 0.0;
	auto arithmetic_mean_L2_error = 0.0;
	auto arithmetic_mean_Linf_error = 0.0;

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		arithmetic_mean_L1_error += this->calculate_cell_L1_error(cell_index, exact_solution, grid, end_time);
		arithmetic_mean_L2_error += this->calculate_cell_L2_error(cell_index, exact_solution, grid, end_time);
		arithmetic_mean_Linf_error += this->calculate_cell_Linf_error(cell_index, exact_solution, grid, end_time);
	}

	arithmetic_mean_L1_error /= num_cells;
	arithmetic_mean_L2_error /= num_cells;
	arithmetic_mean_Linf_error /= num_cells;

	return { arithmetic_mean_L1_error, arithmetic_mean_L2_error, arithmetic_mean_Linf_error };
}

double Semi_Discrete_Equation_DG::calculate_cell_L1_error(const uint cell_index, const Exact_Solution& exact_solution, const Grid& grid, const double end_time) const
{
	const auto solution_degree = this->discrete_solution_->solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = exact_solution.calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = this->discrete_solution_->calculate_solution_at_cell_QPs(cell_index);

	double volume = 0.0;
	double integral_result = 0.0;
	for (ushort q = 0; q < num_QPs; ++q)
	{
		const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
		integral_result += std::abs(error_value) * QWs[q];
		volume += QWs[q];
	}

	return integral_result / volume;
}

double Semi_Discrete_Equation_DG::calculate_cell_L2_error(const uint cell_index, const Exact_Solution& exact_solution, const Grid& grid, const double end_time) const
{
	const auto solution_degree = this->discrete_solution_->solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = exact_solution.calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = this->discrete_solution_->calculate_solution_at_cell_QPs(cell_index);

	double volume = 0.0;
	double integral_result = 0.0;
	for (ushort q = 0; q < num_QPs; ++q)
	{
		const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
		integral_result += error_value * error_value * QWs[q];
		volume += QWs[q];
	}

	return std::sqrt(integral_result / volume);
}

double Semi_Discrete_Equation_DG::calculate_cell_Linf_error(const uint cell_index, const Exact_Solution& exact_solution, const Grid& grid, const double end_time) const
{
	const auto solution_degree = this->discrete_solution_->solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = exact_solution.calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = this->discrete_solution_->calculate_solution_at_cell_QPs(cell_index);

	double max_error = 0.0;
	for (ushort q = 0; q < num_QPs; ++q)
	{
		const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
		max_error = (std::max)(max_error, error_value);
	}

	return max_error;
}

std::unique_ptr<Semi_Discrete_Equation> Semi_Discrete_Equation_Factory::make_unique(const Configuration& configuration, const Grid& grid)
{
	const auto name = configuration.get("spatial_discrete_scheme");

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