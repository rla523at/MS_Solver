#include "../INC/Post_Variables.h"

void Post_Variables::record_grid_data(const Grid& grid, const ushort post_order) {
	this->value_format_convertor_->set_grid_data(grid, post_order);

	const auto set_of_post_nodes = grid.cell_set_of_post_nodes(post_order);
	this->num_post_nodes_ = ms::size_of_vvec(set_of_post_nodes);

	this->set_of_connectivities_ = grid.cell_set_of_connectivities(post_order, set_of_post_nodes);
	this->num_post_elements_ = this->set_of_connectivities_.size();

	this->post_nodes_by_axis_ = this->make_post_coordinate_blocks(set_of_post_nodes);
	this->num_grid_variables_ = this->post_nodes_by_axis_.size();
	REQUIRE(this->num_grid_variables_ <= 3, "num grid variable can not exceed 3");

	if (this->num_grid_variables_ == 2)
		this->zone_type_ = Zone_Type::FETriangle;
	else
		this->zone_type_ = Zone_Type::FETetrahedron;
}

void Post_Variables::record_variable(const std::string_view name, const std::vector<double>& values) {
	REQUIRE(!name.empty(), "post variable should have name");
	REQUIRE(!this->solution_variable_name_to_value_.contains(name.data()), "post variable does not allow duplicate record");

	auto post_variable_values = this->value_format_convertor_->convert_to_post_variable_format(values);
	this->solution_variable_name_to_value_.emplace(name.data(), std::move(post_variable_values));
};

void Post_Variables::set_post_variable_format(const std::string& post_variable_location_str) {
	this->value_format_convertor_ = Post_Variable_Format_Converter_Factory::make(post_variable_location_str);
}

void Post_Variables::syncronize_solution_time(const double& get_solution_time) {
	this->solution_time_ptr_ = &get_solution_time;
}

Zone_Type Post_Variables::zone_type(void) const {
	return this->zone_type_;
}

size_t Post_Variables::num_post_node(void) const {
	return this->num_post_nodes_;
}

size_t Post_Variables::num_post_element(void) const {
	return this->num_post_elements_;
}

ushort Post_Variables::num_grid_variable(void) const {
	return num_grid_variables_;
}

ushort Post_Variables::num_solution_variable(void) const {
	return this->solution_variable_name_to_value_.size();
}

double Post_Variables::solution_time(void) const {
	REQUIRE(this->solution_time_ptr_ != nullptr, "solution time should be syncronized");
	return *solution_time_ptr_;
}

const std::vector<std::vector<double>>& Post_Variables::get_post_nodes_by_axis(void) const {
	return this->post_nodes_by_axis_;
}

const std::vector<std::vector<int>>& Post_Variables::get_connectivities(void) const {
	return this->set_of_connectivities_;
}

std::string Post_Variables::grid_variable_str(void) const {
	if (this->num_grid_variables_ == 1)
		return "x";
	else if (this->num_grid_variables_ == 2)
		return "x,y";
	else if (this->num_grid_variables_ == 3)
		return "x,y,z";
	else
		EXCEPTION("current num grid variable is not supproted");
}

std::vector<std::vector<double>> Post_Variables::make_post_coordinate_blocks(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const {
	const auto space_dimension = set_of_post_nodes.front().front().size();

	std::vector<std::vector<double>> coordinates(space_dimension);

	for (const auto& post_nodes : set_of_post_nodes) {
		for (const auto& node : post_nodes) {
			for (ushort j = 0; j < space_dimension; ++j)
				coordinates[j].push_back(node[j]);
		}
	}

	return coordinates;
}


void Cell_Center_Format_Convertor::set_grid_data(const Grid& grid, const ushort post_order) {
	this->num_post_elements_ = grid.num_cells();

	const auto set_of_num_post_nodes = grid.cell_set_of_num_post_nodes(post_order);
	for (const auto num_post_nodes : set_of_num_post_nodes)
		this->num_post_nodes_ += num_post_nodes;
}

std::vector<double> Cell_Center_Format_Convertor::convert_to_post_variable_format(const std::vector<double>& values) const {
	REQUIRE(values.size() == this->num_post_elements_, "number of values should be same with number of cells");
	return values;
}

std::string Cell_Center_Format_Convertor::solution_variable_location_str(const size_t get_num_solution_variable) const {
	if (get_num_solution_variable == 1)
		return "([1]=CELLCENTERED)";
	else
		return "([1-" + std::to_string(get_num_solution_variable) + "]=CELLCENTERED)";
}

void Node_Format_Convertor::set_grid_data(const Grid& grid, const ushort post_order) {
	this->num_post_elements_ = grid.num_cells();
	this->cell_set_of_num_post_nodes_ = grid.cell_set_of_num_post_nodes(post_order);

	for (const auto num_post_nodes : this->cell_set_of_num_post_nodes_)
		this->num_post_nodes_ += num_post_nodes;
}

std::vector<double> Node_Format_Convertor::convert_to_post_variable_format(const std::vector<double>& values) const {
	const auto num_values = values.size();

	if (this->num_post_nodes_ == num_values)
		return values;
	else {
		REQUIRE(num_values == this->num_post_elements_, "number of values should be same with number of elements");

		std::vector<double> post_variable_values(this->num_post_nodes_);

		size_t index = 0;
		for (size_t i = 0; i < this->num_post_elements_; ++i) {
			const auto num_post_elements = this->cell_set_of_num_post_nodes_[i];
			for (size_t j = 0; j < num_post_elements; ++j)
				post_variable_values[index++] = values[i];
		}

		return post_variable_values;
	}
}

std::string Node_Format_Convertor::solution_variable_location_str(const size_t get_num_solution_variable) const {
	return "()";
}


