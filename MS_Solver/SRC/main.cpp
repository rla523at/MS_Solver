#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Post_Solution_Data.h"
#include "../INC/Log.h"
#include "../INC/Setting.h"

using Grid_Builder_				= Grid_Builder<__DIMENSION__>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX_FUNCTION>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;

int main(void) {
	std::vector<std::string> grid_file_names = { "Quad10","Quad15","Quad20","Quad25","Quad30" };

	for (const auto& grid_file_name : grid_file_names) {
		Log::set_path(__DEFAULT_PATH__ + grid_file_name + "/");
		Post_Solution_Data::set_path(__DEFAULT_PATH__ + grid_file_name + "/");
		Post_Solution_Data::initialize<GOVERNING_EQUATION>(__POST_ORDER__);
		Post_AI_Data::set_path(__DEFAULT_PATH__ + grid_file_name + "/" + "AI_Data/");

		Log::content_ << "================================================================================\n";
		Log::content_ << "\t\t\t\t SETTING \n";
		Log::content_ << "================================================================================\n";
		Log::content_ << std::left << std::setw(35) << "Grid " << grid_file_name << "\n";
		Log::content_ << std::left << std::setw(35) << "Governing Equation" << GOVERNING_EQUATION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Initial Condtion" << INITIAL_CONDITION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Spatial Discrete Method" << SPATIAL_DISCRETE_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Reconstruction Method" << RECONSTRUCTION_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Numeraical Flux Function" << NUMERICAL_FLUX_FUNCTION::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Time Integral Method" << TIME_INTEGRAL_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Time Step Method" << TIME_STEP_METHOD::name() << "\n";
		Log::content_ << std::left << std::setw(35) << "Solve End Condtion" << SOLVE_END_CONDITION::name() << "\n";
		Log::content_ << "================================================================================\n";
		Log::content_ << "================================================================================\n\n";
		Log::print();

		auto grid = Grid_Builder_::build<GRID_FILE_TYPE>(grid_file_name);

		Post_AI_Data::intialize(grid);
		Post_Solution_Data::post_grid(grid.elements.cell_elements);

		Semi_Discrete_Equation_ semi_discrete_equation(std::move(grid));

		auto solutions = semi_discrete_equation.calculate_initial_solutions<INITIAL_CONDITION>();
		Discrete_Equation_::solve<TIME_STEP_METHOD, SOLVE_END_CONDITION, SOLVE_POST_CONDITION>(semi_discrete_equation, solutions);

#ifdef ERROR_CALCULATION
		semi_discrete_equation.estimate_error<INITIAL_CONDITION>(solutions, __END_CONDITION_CONSTANT__);
#endif

		Log::write();
	}
}

