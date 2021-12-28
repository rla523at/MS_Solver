#include "../INC/Error.h"

std::vector<double> No_Error_Calculation::calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const 
{
	return {};
};

Global_Error::Global_Error(std::unique_ptr<Exact_Solution>&& exact_solution)
{
	this->exact_solution_ = std::move(exact_solution);
}

std::vector<double> Global_Error::calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const
{
	const auto L1_norm = this->calculate_global_error_L1_norm(grid, discrete_soltuion, end_time);
	const auto L2_norm = this->calculate_global_error_L2_norm(grid, discrete_soltuion, end_time);
	const auto Linf_norm = this->calculate_global_error_Linf_norm(grid, discrete_soltuion, end_time);

	return { L1_norm, L2_norm, Linf_norm };
}

double Global_Error::calculate_global_error_L1_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto num_cells = grid.num_cells();

	auto total_integral_value = 0.0;
	auto total_volume = 0.0;

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		const auto solution_degree = discrete_solution.solution_degree(cell_index);
		const auto integrand_degree = 2 * solution_degree;
		const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

		const auto& QPs = quadrature_rule.points;
		const auto& QWs = quadrature_rule.weights;
		const auto num_QPs = QPs.size();

		const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
		const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

		auto local_integral_value = 0.0;
		auto local_volume = 0.0;

		for (ushort q = 0; q < num_QPs; ++q)
		{
			const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
			local_integral_value += std::abs(error_value) * QWs[q];
			//local_integral_value += error_value * QWs[q];

			local_volume += QWs[q];
		}

		total_integral_value += local_integral_value;
		//total_integral_value += std::abs(local_integral_value);
		total_volume += local_volume;
	}

	return total_integral_value / total_volume;
}

double Global_Error::calculate_global_error_L2_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto num_cells = grid.num_cells();

	auto total_integral_value = 0.0;
	auto total_volume = 0.0;

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		const auto solution_degree = discrete_solution.solution_degree(cell_index);
		const auto integrand_degree = 2 * solution_degree;
		const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

		const auto& QPs = quadrature_rule.points;
		const auto& QWs = quadrature_rule.weights;
		const auto num_QPs = QPs.size();

		const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
		const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

		auto local_integral_value = 0.0;
		auto local_volume = 0.0;

		for (ushort q = 0; q < num_QPs; ++q)
		{
			const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
			local_integral_value += error_value * error_value * QWs[q];
			local_volume += QWs[q];
		}

		total_integral_value += local_integral_value;
		total_volume += local_volume;
	}

	return std::sqrt(total_integral_value / total_volume);
}

double Global_Error::calculate_global_error_Linf_norm(const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto num_cells = grid.num_cells();

	auto max_error = 0.0;

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		const auto solution_degree = discrete_solution.solution_degree(cell_index);
		const auto integrand_degree = 2 * solution_degree;
		const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

		const auto& QPs = quadrature_rule.points;
		const auto num_QPs = QPs.size();

		const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
		const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

		for (ushort q = 0; q < num_QPs; ++q)
		{
			const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
			max_error = (std::max)(max_error, std::abs(error_value));
		}
	}

	return max_error;
}

Local_Average_Error::Local_Average_Error(std::unique_ptr<Exact_Solution>&& exact_solution)
{
	this->exact_solution_ = std::move(exact_solution);
}

std::vector<double> Local_Average_Error::calculate_error_norms(const Grid& grid, const Discrete_Solution_DG& discrete_soltuion, const double end_time) const
{
	const auto num_cells = grid.num_cells();

	auto arithmetic_mean_L1_error = 0.0;
	auto arithmetic_mean_L2_error = 0.0;
	auto arithmetic_mean_Linf_error = 0.0;

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		arithmetic_mean_L1_error += this->calculate_local_error_L1_norm(cell_index, grid, discrete_soltuion, end_time);
		arithmetic_mean_L2_error += this->calculate_local_error_L2_norm(cell_index, grid, discrete_soltuion, end_time);
		arithmetic_mean_Linf_error += this->calculate_local_error_Linf_norm(cell_index, grid, discrete_soltuion, end_time);
	}

	arithmetic_mean_L1_error /= num_cells;
	arithmetic_mean_L2_error /= num_cells;
	arithmetic_mean_Linf_error /= num_cells;

	return { arithmetic_mean_L1_error, arithmetic_mean_L2_error, arithmetic_mean_Linf_error };
}

double Local_Average_Error::calculate_local_error_L1_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto solution_degree = discrete_solution.solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

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

double Local_Average_Error::calculate_local_error_L2_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto solution_degree = discrete_solution.solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

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

double Local_Average_Error::calculate_local_error_Linf_norm(const uint cell_index, const Grid& grid, const Discrete_Solution_DG& discrete_solution, const double end_time) const
{
	const auto solution_degree = discrete_solution.solution_degree(cell_index);
	const auto integrand_degree = 2 * solution_degree;
	const auto& quadrature_rule = grid.get_cell_quadrature_rule(cell_index, integrand_degree);

	const auto& QPs = quadrature_rule.points;
	const auto& QWs = quadrature_rule.weights;
	const auto num_QPs = QPs.size();

	const auto exact_solution_v_at_QPs = this->exact_solution_->calculate_exact_solution_vectors(QPs, end_time);
	const auto computed_solution_v_at_QPs = discrete_solution.calculate_solution_at_cell_QPs(cell_index);

	double max_error = 0.0;
	for (ushort q = 0; q < num_QPs; ++q)
	{
		const auto error_value = (exact_solution_v_at_QPs[q] - computed_solution_v_at_QPs[q])[0];
		max_error = (std::max)(max_error, error_value);
	}

	return max_error;
}

std::unique_ptr<Error> Error_Factory::make_unqiue(const Configuration& configuration)
{
	const auto& write_error_file = configuration.get_write_error_file();

	if (ms::compare_icase(write_error_file, "No"))
	{
		return std::make_unique<No_Error_Calculation>();
	}

	const auto& governing_equation_name = configuration.get_governing_equation();

	if (ms::contains_icase(governing_equation_name, "Linear", "Advection"))
	{
		const auto& error_type = configuration.get_error_type();
		auto exact_solution = Exact_Solution_Factory::make_unique(configuration);

		if (ms::contains_icase(error_type, "Global"))
		{
			return std::make_unique<Global_Error>(std::move(exact_solution));
		}
		else if (ms::contains_icase(error_type, "Local_Average"))
		{
			return std::make_unique<Local_Average_Error>(std::move(exact_solution));
		}
		else
		{
			EXCEPTION("error type in configuration file is not supported");
			return nullptr;
		}
	}
	else
	{
		return std::make_unique<No_Error_Calculation>();
	}
}