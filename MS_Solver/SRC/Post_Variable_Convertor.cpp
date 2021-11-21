#include "../INC/Post_Variable_Convertor.h"

Post_Variable_Convertor::Post_Variable_Convertor(const Grid& grid, const ushort post_order, const Discrete_Solution& discrete_solution)
	:discrete_solution_(discrete_solution)
{
	this->num_cells_ = grid.num_cells();
}

const std::vector<std::string>& Post_Variable_Convertor::get_solution_names(void) const
{
	return this->discrete_solution_.get_solution_names();
}

Node_Base_Convertor::Node_Base_Convertor(const Grid& grid, const ushort post_order, Discrete_Solution& discrete_solution)
	: Post_Variable_Convertor(grid, post_order, discrete_solution) 
{
	this->set_of_num_post_points_ = grid.cell_set_of_num_post_points(post_order);

	for (int i = 0; i < this->num_cells_; ++i)
	{
		this->num_post_points_ += this->set_of_num_post_points_[i];
	}

	discrete_solution.precalculate_post_points(grid.cell_set_of_post_points(post_order));
};

std::vector<double> Node_Base_Convertor::convert_values(const std::vector<double>& values) const
{
	REQUIRE(values.size() == this->num_cells_, "number of values should be same with number of cells");

	std::vector<double> post_variable_values(this->num_post_points_);

	size_t index = 0;
	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto num_post_elements = this->set_of_num_post_points_[cell_index];

		for (int j = 0; j < num_post_elements; ++j)
		{
			post_variable_values[index++] = values[cell_index];
		}
	}

	return post_variable_values;
}

std::vector<std::vector<double>> Node_Base_Convertor::calculate_set_of_post_point_solution_values(void) const
{
	const auto num_solutions = this->discrete_solution_.num_solutions();

	std::vector<std::vector<double>> set_of_post_point_solution_values(num_solutions);

	for (auto& solution_values : set_of_post_point_solution_values)
	{
		solution_values.reserve(this->num_post_points_);
	}

	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto solution_at_post_points = this->discrete_solution_.calculate_solution_at_post_points(cell_index);

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

std::string Node_Base_Convertor::variable_location_str(void) const 
{
	return "Node";
}

Center_Base_Convertor::Center_Base_Convertor(const Grid& grid, const ushort post_order, Discrete_Solution& discrete_solution)
	: Post_Variable_Convertor(grid, post_order, discrete_solution)
{
	this->set_of_num_post_elements_ = grid.cell_set_of_num_post_elements(post_order);

	for (int i = 0; i < this->num_cells_; ++i)
	{
		this->num_post_elements_ += this->set_of_num_post_elements_[i];
	}

	discrete_solution.precalculate_post_elements(grid.cell_set_of_post_element_centers(post_order));
};

std::vector<double> Center_Base_Convertor::convert_values(const std::vector<double>& values) const
{
	REQUIRE(values.size() == this->num_cells_, "number of values should be same with number of cells");

	std::vector<double> post_variable_values(this->num_post_elements_);

	size_t index = 0;
	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto num_post_elements = this->set_of_num_post_elements_[cell_index];
		for (int j = 0; j < num_post_elements; ++j)
		{
			post_variable_values[index++] = values[cell_index];
		}
	}

	return post_variable_values;
}

std::vector<std::vector<double>> Center_Base_Convertor::calculate_set_of_post_point_solution_values(void) const
{
	const auto num_solutions = this->discrete_solution_.num_solutions();

	std::vector<std::vector<double>> set_of_post_point_solution_values(num_solutions);

	for (auto& solution_values : set_of_post_point_solution_values)
	{
		solution_values.reserve(this->num_post_elements_);
	}

	for (int cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto solution_at_post_element_centers = this->discrete_solution_.calculate_solution_at_post_element_centers(cell_index);

		for (int i = 0; i < this->set_of_num_post_elements_[cell_index]; ++i)
		{
			const auto& solution = solution_at_post_element_centers[i];

			for (int j = 0; j < num_solutions; ++j)
			{
				set_of_post_point_solution_values[j].push_back(solution[j]);
			}
		}
	}

	return set_of_post_point_solution_values;
}

std::string Center_Base_Convertor::variable_location_str(void) const
{
	return "Center";
}

std::unique_ptr<Post_Variable_Convertor> Post_Variable_Converter_Factory::make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution& discrete_solution)
{
	const auto post_order = configuration.get<ushort>("post_order");
	const auto post_point_location = configuration.get("post_point_location");

	if (ms::contains_icase(post_point_location, "Node"))
	{
		return std::make_unique<Node_Base_Convertor>(grid, post_order, discrete_solution);
	}
	else if (ms::contains_icase(post_point_location, "Center"))
	{
		return std::make_unique<Center_Base_Convertor>(grid, post_order, discrete_solution);
	}
	else
	{
		EXCEPTION("post point loacation in configuration file does not supported");
	}
}
