#include "../INC/Inital_Condition.h"
#include "../INC/Discrete_Equation.h"
#include "../INC/Post_Solution_Data.h"
#include "../INC/Log.h"
#include "../INC/Setting.h"

using Post_Solution_Data_		= Post_Solution_Data<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, __POST_ORDER__>;
using Grid_Builder_				= Grid_Builder<__DIMENSION__>;
using Semi_Discrete_Equation_	= Semi_Discrete_Equation<GOVERNING_EQUATION, SPATIAL_DISCRETE_METHOD, RECONSTRUCTION_METHOD, NUMERICAL_FLUX_FUNCTION>;
using Discrete_Equation_		= Discrete_Equation<TIME_INTEGRAL_METHOD>;

int main(void) {
	Log::set_path(PATH);
	Post_Solution_Data_::set_path(PATH);
	Post_AI_Data::set_path(PATH + "AI_Data/");

	Log::content_ << "================================================================================\n";
	Log::content_ << "\t\t\t\t SETTING \n";
	Log::content_ << "================================================================================\n";
	Log::content_ << std::left << std::setw(35) << "Grid " << GRID_FILE_NAME << "\n";
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

	auto grid = Grid_Builder_::build<GRID_FILE_TYPE>(GRID_FILE_NAME);
		
	Post_AI_Data::intialize(grid);
	Post_Solution_Data_::post_grid(grid.elements.cell_elements);

	auto semi_discrete_equation = Semi_Discrete_Equation_(std::move(grid));

	auto solutions = semi_discrete_equation.calculate_initial_solutions<INITIAL_CONDITION>();
	Discrete_Equation_::solve<TIME_STEP_METHOD, SOLVE_END_CONDITION, SOLVE_POST_CONDITION, Post_Solution_Data_>(semi_discrete_equation, solutions);
		
#ifdef ERROR_CALCULATION
	Semi_Discrete_Equation_::estimate_error<INITIAL_CONDITION>(solutions, __END_CONDITION_CONSTANT__);
#endif

	Log::write();
}