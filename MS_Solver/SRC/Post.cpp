#include "../INC/Post.h"

void Post_Processing::initialize(const Configuration& configuration) {	
	This_::post_file_path_ = configuration.get("Post_file_path");

	This_::post_order_ = configuration.get<ushort>("post_order");

	const auto post_variable_location = configuration.get("post_variable_location");	
	This_::post_variables_.set_post_variable_format(post_variable_location);

	if (ms::contains_icase(post_variable_location, "cell", "center"))
		This_::post_order_ = 0;
};

void Post_Processing::post_grid(const Grid& grid) {
	This_::post_variables_.record_grid_data(grid, This_::post_order_);
	This_::file_writer_.write_grid_file(This_::post_variables_, This_::post_file_path_);
}

void Post_Processing::post_solution(const Discretized_Solution& discretized_solution) {
	const auto& variable_names = discretized_solution.get_variable_names();
	const auto post_point_solutions_by_variable = discretized_solution.calculate_post_point_solutions_by_variable();
	const auto get_num_solution_variable = variable_names.size();

	for (ushort i = 0; i < get_num_solution_variable; ++i)
		This_::post_variables_.record_variable(variable_names[i], post_point_solutions_by_variable[i]);

	This_::file_writer_.write_solution_file(This_::post_variables_, post_file_path_);
}

void Post_Processing::record_variables(const std::string_view name, const std::vector<double>& values) {
	This_::post_variables_.record_variable(name, values);
}

void Post_Processing::syncronize_solution_time(const double& get_solution_time) {
	This_::post_variables_.syncronize_solution_time(get_solution_time);
}
