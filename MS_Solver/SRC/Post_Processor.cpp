#include "../INC/Post_Processor.h"

void Post_Processor::initialize(const Configuration& configuration, const Grid& grid, const Discrete_Solution_DG& discrete_solution) 
{	
	This_::post_file_path_ = configuration.get("Post_file_path");
	This_::post_order_ = configuration.get<ushort>("post_order");

	auto post_variable_convertor = Post_Variable_Converter_Factory::make_unique(grid, This_::post_order_, discrete_solution);
	This_::post_variables_ = std::make_unique<Post_Variables>(grid, std::move(post_variable_convertor), This_::post_order_);
	This_::file_writer_ = Tecplot_File_Writer_Factory::make_unique(configuration);

	This_::is_initialized_ = true;
};

void Post_Processor::post_grid(void) 
{
	REQUIRE(This_::is_initialized_, "Post processore should be initialized before post grid");

	This_::file_writer_->write_grid_file(*This_::post_variables_, This_::post_file_path_);
}

void Post_Processor::post_solution(void) 
{
	REQUIRE(This_::is_initialized_, "Post processore should be initialized before post solution");
	
	This_::post_variables_->record_solution();
	This_::file_writer_->write_solution_file(*This_::post_variables_, post_file_path_);
}

void Post_Processor::record_variables(const std::string& name, const std::vector<double>& values) 
{
	REQUIRE(This_::is_initialized_, "Post processore should be initialized before record variable");

	This_::post_variables_->record_variable(name, values);
}

void Post_Processor::syncronize_solution_time(const double& solution_time) 
{
	REQUIRE(This_::is_initialized_, "Post processore should be initialized before syncronize solution time");

	This_::post_variables_->syncronize_solution_time(solution_time);
}
