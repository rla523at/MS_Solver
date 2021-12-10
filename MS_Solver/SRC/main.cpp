#include "../INC/Configuration.h"
#include "../INC/Discrete_Equation.h"

int main(void) 
{	
	try
	{
		const auto configuration_file_path = "RSC/configuration.dat";
		Configuration configuration(configuration_file_path);
		
		const auto& grid_file_paths = configuration.get_grid_file_paths();
		for (const auto& grid_file_path : grid_file_paths)
		{
			const auto grid_file_name = configuration.grid_file_name(grid_file_path);
			LOG << configuration.configuration_str(grid_file_name);

			const auto post_folder_path = configuration.post_folder_path_str(grid_file_name);
			Post_Processor::set_path(post_folder_path);
			LOG.set_path(post_folder_path);			

			LOG << "================================================================================\n";
			LOG << "\t\t\t\t Start Pre-Processing \n";
			LOG << "================================================================================\n" << Log::print_;

			Profiler::set_time_point();

			auto grid_file_convertor = Grid_File_Convertor_Factory::make_unique(configuration);
			Grid grid(*grid_file_convertor, grid_file_path);

			auto semi_discrete_equation = Semi_Discrete_Equation_Factory::make_unique(configuration, grid);
			auto time_discrete_scheme = Time_Discrete_Scheme_Factory::make_unique(configuration);
			auto solve_end_controller = Solve_End_Controller_Factory::make_unique(configuration);
			auto solve_post_controller = Solve_Post_Controller_Factory::make_unique(configuration);

			Discrete_Equation discrete_equation(std::move(semi_discrete_equation), std::move(time_discrete_scheme), std::move(solve_end_controller), std::move(solve_post_controller));

			LOG << "\n================================================================================\n";
			LOG << "\t\t\t End Pre-Processing(" << Profiler::get_time_duration() << "s)\n";
			LOG << "================================================================================\n\n" << LOG.print_;

			LOG << "================================================================================\n";
			LOG << "\t\t\t\t Start Solving\n";
			LOG << "================================================================================\n" << LOG.print_;

			Profiler::set_time_point();

			discrete_equation.solve();

			LOG << "\n================================================================================\n";
			LOG << "\t\t\t End Solving(" << Profiler::get_time_duration() << "s)\n";
			LOG << "================================================================================\n\n" << LOG.print_;

			LOG << "================================================================================\n";
			LOG << "\t\t\t\t Start Error Analysis\n";
			LOG << "================================================================================\n";

			const auto exact_solution = Exact_Solution_Factory::make_unique(configuration);

			Profiler::set_time_point();

			//const auto error_values = discrete_equation.calculate_global_error_values(*exact_solution, grid);
			const auto error_values = discrete_equation.calculate_specific_error_values(*exact_solution, grid);

			LOG << std::left << std::setprecision(16);
			LOG << std::setw(25) << "L1 error" << std::setw(25) << "L2 error" << "Linf error \n";
			LOG << std::setw(25) << error_values[0] << std::setw(25) << error_values[1] << error_values[2] << "\n";

			LOG << "\n================================================================================\n";
			LOG << "\t\t\t End Error Analysis(" << Profiler::get_time_duration() << "s)\n";
			LOG << "================================================================================\n\n" << LOG.print_;

			LOG.write_error_text(grid_file_name, error_values);
			LOG.write();
			LOG.clear();
		}
	}
	catch (const std::exception& except)
	{
		LOG << "\n================================================================================\n";
		LOG << except.what() << "\n";
		LOG << "================================================================================\n\n" << LOG.print_;
	}
}