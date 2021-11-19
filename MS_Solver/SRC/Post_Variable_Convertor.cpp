#include "../INC/Post_Variable_Convertor.h"

Post_Variable_Convertor::Post_Variable_Convertor(const Grid& grid, const ushort post_order)
{
	this->num_cells_ = grid.num_cells();
	this->set_of_num_post_elements_ = grid.cell_set_of_num_post_elements(post_order);
	this->set_of_num_post_points_ = grid.cell_set_of_num_post_points(post_order);

	for (int i = 0; i < this->num_cells_; ++i)
	{
		this->num_post_elements_ += this->set_of_num_post_elements_[i];
		this->num_post_points_ += this->set_of_num_post_points_[i];
	}
}

std::vector<double> Post_Variable_Convertor::convert_cell_center_values(const std::vector<double>& values) const
{
	REQUIRE(values.size() == this->num_cells_, "number of values should be same with number of cells");

	std::vector<double> post_variable_values(this->num_post_elements_);

	size_t index = 0;
	for (size_t i = 0; i < this->num_cells_; ++i)
	{
		const auto num_post_elements = this->set_of_num_post_elements_[i];
		for (size_t j = 0; j < num_post_elements; ++j)
		{
			post_variable_values[index++] = values[i];
		}
	}

	return post_variable_values;
}

std::vector<double> Post_Variable_Convertor::convert_values(const std::vector<double>& values) const
{
	REQUIRE(values.size() == this->num_cells_, "number of values should be same with number of cells");

	std::vector<double> post_variable_values(this->num_post_points_);

	size_t index = 0;
	for (size_t i = 0; i < this->num_cells_; ++i)
	{
		const auto num_post_elements = this->set_of_num_post_points_[i];
		for (size_t j = 0; j < num_post_elements; ++j)
		{
			post_variable_values[index++] = values[i];
		}
	}

	return post_variable_values;
}

Post_Variable_Convertor_DG::Post_Variable_Convertor_DG(const Grid& grid, const ushort post_order, const Discrete_Solution_DG& discrete_solution)
	: Post_Variable_Convertor(grid, post_order), discrete_solution_(discrete_solution)
{
	this->P0_basis_values_.resize(this->num_cells_);
	this->set_of_basis_post_points_m_.reserve(this->num_cells_);

	const auto set_of_post_points = grid.cell_set_of_post_points(post_order);
	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		this->P0_basis_values_[cell_index] = this->discrete_solution_.calculate_P0_basis_value(cell_index);

		const auto& post_points = set_of_post_points[cell_index];
		auto basis_post_points_m = this->discrete_solution_.calculate_basis_points_m(cell_index, post_points);

		this->set_of_basis_post_points_m_.push_back(std::move(basis_post_points_m));
	}
}

std::vector<std::vector<double>> Post_Variable_Convertor_DG::calculate_set_of_post_point_solution_values(void) const
{
	const auto num_solutions = this->discrete_solution_.num_solutions();

	std::vector<std::vector<double>> set_of_post_point_solution_values(num_solutions);

	for (auto& solution_values : set_of_post_point_solution_values)
	{
		solution_values.reserve(this->num_post_points_);
	}

	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto solution_at_post_points = this->discrete_solution_.calculate_solution_at_points(cell_index, this->set_of_basis_post_points_m_[cell_index]);

		for (int j = 0; j < this->set_of_num_post_points_[cell_index]; ++j)
		{
			const auto& solution = solution_at_post_points[j];

			for (int i = 0; i < num_solutions; ++i)
			{
				set_of_post_point_solution_values[i].push_back(solution[i]);
			}
		}
	}

	return set_of_post_point_solution_values;
}

std::vector<std::vector<double>> Post_Variable_Convertor_DG::calculate_set_of_cell_center_solution_values(void) const
{
	const auto num_solutions = this->discrete_solution_.num_solutions();

	std::vector<std::vector<double>> set_of_post_point_solution_values(num_solutions);

	for (auto& solution_values : set_of_post_point_solution_values)
	{
		solution_values.reserve(this->num_post_elements_);
	}

	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto P0_solution = this->discrete_solution_.calculate_P0_solution(cell_index, this->P0_basis_values_[cell_index]);

		for (int i = 0; i < this->set_of_num_post_elements_[cell_index]; ++i)
		{
			for (int j = 0; j < num_solutions; ++j)
			{
				set_of_post_point_solution_values[j].push_back(P0_solution[j]);
			}
		}
	}

	return set_of_post_point_solution_values;
}

const std::vector<std::string>& Post_Variable_Convertor_DG::get_solution_names(void) const
{
	return this->discrete_solution_.get_solution_names();
}

std::unique_ptr<Post_Variable_Convertor> Post_Variable_Converter_Factory::make_unique(const Grid& grid, const ushort post_order, const Discrete_Solution_DG& discrete_solution)
{
	return std::make_unique<Post_Variable_Convertor_DG>(grid, post_order, discrete_solution);
}