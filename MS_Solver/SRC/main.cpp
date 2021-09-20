#include "../INC/Setting.h"

using Grid_Element_Builder_		= Grid_Element_Builder<GRID_FILE_TYPE, __DIMENSION__>;
using Grid_						= Grid<__DIMENSION__>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX_FUNCTION, SCAILING_METHOD_FLAG>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;

int main(void) {
	ms::apply_user_defined_setting();

	const char delimiter = ',';
	std::string grid_file_names_str = GRID_FILE_NAMES;
	auto grid_file_names = ms::parse(grid_file_names_str, delimiter);
	

	for (auto& grid_file_name : grid_file_names) {
		ms::be_removed(grid_file_name, " ");

		Tecplot::initialize<GOVERNING_EQUATION>(__POST_ORDER__, POST_FILE_FORMAT); //post
		Solve_Controller::initialize(SOLVE_END_CONDITION, SOLVE_POST_CONDITION);

		const auto date_str = Log::date_string();

		Log::content_ << "current date : " << date_str << "\n\n";
		Log::content_ << "================================================================================\n";
		Log::content_ << "\t\t\t\t SETTING \n";
		Log::content_ << "================================================================================\n";
		Log::content_ << std::left << std::setw(35) << "Grid " << grid_file_name << "\n";
		Log::content_ << std::left << std::setw(35) << "Governing Equation" << GOVERNING_EQUATION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Initial Condtion" << INITIAL_CONDITION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Spatial Discrete Method" << SPATIAL_DISCRETE_METHOD::name() << "\n";

		if constexpr (SCAILING_METHOD_FLAG)
			Log::content_ << std::left << std::setw(35) << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << " with scailing method\n";
		else
			Log::content_ << std::left << std::setw(35) << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << "\n";

		Log::content_ << std::left << std::setw(35) << "Numeraical Flux Function" << NUMERICAL_FLUX_FUNCTION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Time Integral Method" << TIME_INTEGRAL_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Time Step Method" << TIME_STEP_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Solve End Condtion" << Solve_Controller::end_condition_name() << "\n";
		Log::content_ << "================================================================================\n";
		Log::content_ << "================================================================================\n\n";
		Log::print();

		Log::set_path(__DEFAULT_PATH__ + grid_file_name + "_" + date_str + "/");
		Tecplot::set_path(__DEFAULT_PATH__ + grid_file_name + "_" + date_str + "/"); //post

		auto grid_element = Grid_Element_Builder_::build_from_grid_file(grid_file_name);
		Grid_ grid(std::move(grid_element));

		Tecplot::post_grid(grid); //post

		Semi_Discrete_Equation_ semi_discrete_equation(grid);
		auto solutions = semi_discrete_equation.calculate_initial_solutions<INITIAL_CONDITION>();


		try { Discrete_Equation_::solve<TIME_STEP_METHOD>(semi_discrete_equation, solutions); }
		catch (const std::exception& exception) {
			Log::content_ << "\n================================================================================\n";
			Log::content_ << "\t\t\t Abnormal Termination\n";
			Log::content_ << "================================================================================\n";
			Log::content_ << "Essential requirement is not satisfied => " << exception.what() << "\n";
			Log::print();
			Log::write();
			Tecplot::post_solution(solutions, "abnormal_termination");
			std::exit(523);
		}

		semi_discrete_equation.estimate_error<INITIAL_CONDITION>(solutions, __SOLVE_END_CONDITION_CONSTANT__);

		Log::write();
	}
}

