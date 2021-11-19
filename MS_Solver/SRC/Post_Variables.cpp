#include "../INC/Post_Variables.h"

Post_Variables::Post_Variables(const Grid& grid, std::unique_ptr<Post_Variable_Convertor>&& post_variable_convertor, const ushort post_order)
	:post_variable_convertor_(std::move(post_variable_convertor))
{
	this->space_dimension_ = grid.space_dimension();

	const auto set_of_post_points = grid.cell_set_of_post_points(post_order);
	this->num_post_points_ = ms::size_of_vvec(set_of_post_points);

	this->set_of_connectivities_ = grid.cell_set_of_connectivities(post_order, set_of_post_points);
	this->num_post_elements_ = this->set_of_connectivities_.size();

	this->post_coordinate_blocks = this->make_post_coordinate_blocks(set_of_post_points);

	if (this->space_dimension_ == 2)
	{
		this->zone_type_ = Zone_Type::FETriangle;
	}
	else if (this->space_dimension_ == 3)
	{
		this->zone_type_ = Zone_Type::FETetrahedron;
	}
	else
	{
		EXCEPTION("not supported space dimension");
	}
}

void Post_Variables::record_cell_center_solution(void)
{
	const auto& solution_names = this->post_variable_convertor_->get_solution_names();
	auto set_of_cell_center_solution_values = this->post_variable_convertor_->calculate_set_of_cell_center_solution_values();
	const auto num_solutions = solution_names.size();

	for (ushort i = 0; i < num_solutions; ++i)
	{
		this->cell_center_post_variable_names_.push_back(solution_names[i]);
		this->set_of_cell_center_post_variable_values_.push_back(std::move(set_of_cell_center_solution_values[i]));
	}
}

void Post_Variables::record_cell_center_variable(const std::string& name, const std::vector<double>& values)
{
	REQUIRE(!name.empty(), "post variable should have name");
	REQUIRE(!ms::contains(post_variable_names_, name), "post variable does not allow duplicate record");

	auto post_variable_values = this->post_variable_convertor_->convert_cell_center_values(values);

	this->cell_center_post_variable_names_.push_back(name);
	this->set_of_cell_center_post_variable_values_.push_back(std::move(post_variable_values));
}

void Post_Variables::record_solution(void)
{
	const auto& solution_names = this->post_variable_convertor_->get_solution_names();
	auto set_of_post_point_solution_values = this->post_variable_convertor_->calculate_set_of_post_point_solution_values();
	const auto num_solutions = solution_names.size();

	for (ushort i = 0; i < num_solutions; ++i)
	{
		this->post_variable_names_.push_back(solution_names[i]);
		this->set_of_post_variable_values_.push_back(std::move(set_of_post_point_solution_values[i]));
	}
}

void Post_Variables::record_variable(const std::string& name, const std::vector<double>& values)
{
	REQUIRE(!name.empty(), "post variable should have name");
	REQUIRE(!ms::contains(post_variable_names_, name), "post variable does not allow duplicate record");

	auto post_variable_values = this->post_variable_convertor_->convert_values(values);

	this->post_variable_names_.push_back(name);
	this->set_of_post_variable_values_.push_back(std::move(post_variable_values));
}

void Post_Variables::syncronize_solution_time(const double& solution_time) 
{
	this->solution_time_ptr_ = &solution_time;
}

const std::vector<std::vector<double>>& Post_Variables::get_post_coordinate_blocks(void) const
{
	return this->post_coordinate_blocks;
}

const std::vector<std::vector<int>>& Post_Variables::get_connectivities(void) const
{
	return this->set_of_connectivities_;
}

const std::vector<std::vector<double>>& Post_Variables::get_set_of_post_variable_values(void) const
{
	return this->set_of_post_variable_values_;
}

const std::vector<std::string>& Post_Variables::get_post_variable_names(void) const
{
	return this->post_variable_names_;
}

const std::vector<std::vector<double>>& Post_Variables::get_set_of_cell_center_post_variable_values(void) const
{
	return this->set_of_cell_center_post_variable_values_;
}

const std::vector<std::string>& Post_Variables::get_cell_center_post_variable_names(void) const
{
	return this->cell_center_post_variable_names_;
}

std::string Post_Variables::grid_variable_str(void) const
{
	if (this->space_dimension_ == 1)
	{
		return "x";
	}
	else if (this->space_dimension_ == 2)
	{
		return "x,y";
	}
	else if (this->space_dimension_ == 3)
	{
		return "x,y,z";
	}
	else
	{
		EXCEPTION("current num grid variable is not supproted");
		return "";
	}
}

size_t Post_Variables::num_post_node(void) const 
{
	return this->num_post_points_;
}

size_t Post_Variables::num_post_element(void) const 
{
	return this->num_post_elements_;
}

ushort Post_Variables::num_grid_variable(void) const 
{
	return this->space_dimension_;
}

ushort Post_Variables::num_post_variables(void) const 
{
	return this->post_variable_names_.size();
}

ushort Post_Variables::num_cell_center_post_variables(void) const
{
	return this->cell_center_post_variable_names_.size();
}

double Post_Variables::solution_time(void) const 
{
	REQUIRE(this->solution_time_ptr_ != nullptr, "solution time should be syncronized");
	return *solution_time_ptr_;
}

Zone_Type Post_Variables::zone_type(void) const
{
	return this->zone_type_;
}

std::vector<std::vector<double>> Post_Variables::make_post_coordinate_blocks(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points) const 
{
	std::vector<std::vector<double>> coordinate_blocks(this->space_dimension_);
	
	for (auto& coordinate_block : coordinate_blocks)
	{
		coordinate_block.reserve(this->num_post_points_);
	}

	for (const auto& post_points : set_of_post_points)
	{
		for (const auto& point : post_points) 
		{
			for (int i = 0; i < this->space_dimension_; ++i)
			{
				coordinate_blocks[i].push_back(point[i]);
			}
		}
	}

	return coordinate_blocks;
}

