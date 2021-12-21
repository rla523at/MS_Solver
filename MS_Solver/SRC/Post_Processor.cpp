#include "../INC/Post_Processor.h"
void Post_Processor::set_configuration(const Configuration& configuration)
{
	const auto& write_post_file = configuration.get_write_post_file();
	if (ms::compare_icase(write_post_file, "Yes"))
	{
		This_::do_write_ = true;
	}
	else if (ms::compare_icase(write_post_file, "No"))
	{
		This_::do_write_ = false;
	}
	else
	{
		EXCEPTION("post processor switch in configuaration file is not supported");
	}
}

void Post_Processor::initialize(const Configuration& configuration, const Grid& grid, Discrete_Solution& discrete_solution) 
{	
	Profiler::set_time_point();

	This_::post_order_ = configuration.post_order();

	auto post_variable_convertor = Post_Variable_Converter_Factory::make_unique(configuration, grid, discrete_solution);
	This_::post_variables_ = std::make_unique<Post_Variables>(grid, std::move(post_variable_convertor), This_::post_order_);
	This_::file_writer_ = Tecplot_File_Writer_Factory::make_unique(configuration);

	REQUIRE(!This_::post_folder_path_.empty(), "post folder path should be set");
	This_::is_initialized_ = true;

	LOG << std::left << std::setw(50) << "@ Post Processor precalculation" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << LOG.print_;
};

void Post_Processor::post_grid(void) 
{
	if (!This_::do_write_)
	{
		return;
	}

	REQUIRE(This_::is_initialized_, "Post processore should be initialized before post grid");

	const auto post_file_path = This_::post_folder_path_ + "grid.plt";
	This_::file_writer_->write_grid_file(*This_::post_variables_, post_file_path);
}


void Post_Processor::post_solution(const std::string_view sv)
{
	if (!This_::do_write_)
	{
		return;
	}

	REQUIRE(This_::is_initialized_, "Post processore should be initialized before post solution");
	REQUIRE(This_::is_solution_recorded_, "Post processore should be recorded before post solution");

	const auto post_file_name = "solution" + std::to_string(This_::file_writer_->strand_id()) + "_" + sv.data() + ".plt";
	const auto post_file_path = This_::post_folder_path_ + post_file_name;
	This_::file_writer_->write_solution_file(*This_::post_variables_, post_file_path);

	This_::post_variables_->clear_variables();
	This_::is_solution_recorded_ = false;
}

void Post_Processor::record_solution(void)
{
	if (!This_::do_write_)
	{
		return;
	}

	This_::post_variables_->record_solution();
	This_::is_solution_recorded_ = true;
}

void Post_Processor::record_variables(const std::string& name, const std::vector<double>& values) 
{
	if (!This_::do_write_)
	{
		return;
	}

	REQUIRE(This_::is_initialized_, "Post processore should be initialized before record variable");

	This_::post_variables_->record_variable(name, values);
}

void Post_Processor::set_path(const std::string& path)
{
	This_::post_folder_path_ = path;
}

void Post_Processor::syncronize_solution_time(const double& solution_time) 
{
	if (!This_::do_write_)
	{
		return;
	}

	REQUIRE(This_::is_initialized_, "Post processore should be initialized before syncronize solution time");

	This_::post_variables_->syncronize_solution_time(solution_time);
}

