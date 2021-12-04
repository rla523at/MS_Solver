#include "gtest/gtest.h"
#include "../MS_Solver/INC/Discrete_Solution.h"
#include "../MS_Solver/INC/Inner_Faces.h"
#include "../MS_Solver/INC/Post_Processor.h"

TEST(Discrete_Solution_DG, turn_off_log_print)
{
	LOG << Log::print_off;
}
TEST(Discrete_Solution_DG, calculate_solution_at_cell_QPs_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/RQ3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	for (int solution_degree = 0; solution_degree < 8; ++solution_degree)
	{
		const auto governing_equation = std::make_shared<Burgers_2D>();
		const auto initial_condition = std::make_unique<Constant1>();
		auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		const auto num_cells = grid.num_cells();
		std::vector<Quadrature_Rule> cell_quadrature_rules(num_cells);

		for (int cell_index = 0; cell_index < num_cells; ++cell_index)
		{
			cell_quadrature_rules[cell_index] = grid.get_cell_quadrature_rule(cell_index, solution_degree);
		}

		discrete_solution_DG->precalcualte_cell_QPs_basis_values(cell_quadrature_rules);


		for (int cell_index = 0; cell_index < num_cells; ++cell_index)
		{
			const auto solution_at_cell_QPs = discrete_solution_DG->calculate_solution_at_cell_QPs(cell_index);

			for (const auto& solution : solution_at_cell_QPs)
			{
				constexpr auto ref_solution = 1.0;
				constexpr auto epsilon = 9.0E-12;
				EXPECT_NEAR(ref_solution, solution[0], epsilon);
				//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
			}
		}
	}
}
TEST(Discrete_Solution_DG, calculate_P0_solutions_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/RQ3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	for (int solution_degree = 0; solution_degree < 8; ++solution_degree)
	{
		const auto governing_equation = std::make_shared<Burgers_2D>();
		const auto initial_condition = std::make_unique<Constant1>();
		auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		discrete_solution_DG->precalculate_cell_P0_basis_values();

		const auto P0_solutions = discrete_solution_DG->calculate_P0_solutions();
		for (const auto& solution : P0_solutions)
		{
			constexpr auto ref_solution = 1.0;
			constexpr auto epsilon = 9.0E-15;
			EXPECT_NEAR(ref_solution, solution[0], epsilon);
			//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
		}
	}
}
TEST(Discrete_Solution_DG, calculate_solution_at_infc_QPs_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/RQ3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	for (int solution_degree = 0; solution_degree < 8; ++solution_degree)
	{
		const auto governing_equation = std::make_shared<Burgers_2D>();
		const auto initial_condition = std::make_unique<Constant1>();
		auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		const auto numerical_flux_function = std::make_shared<LLF>(governing_equation);
		auto inner_faces_DG = std::make_unique<Inner_Faces_DG>(numerical_flux_function, grid, *discrete_solution_DG);

		const auto num_infcs = grid.num_inner_faces();
		for (int infc_index = 0; infc_index < num_infcs; ++infc_index)
		{
			const auto [oc_index, nc_index] = grid.inner_face_oc_nc_index_pair(infc_index);

			const auto solution_at_ocs_infc_QPs = discrete_solution_DG->calculate_solution_at_infc_ocs_QPs(infc_index, oc_index);
			const auto solution_at_ncs_infc_QPs = discrete_solution_DG->calculate_solution_at_infc_ncs_QPs(infc_index, nc_index);
			
			for (const auto& solution : solution_at_ocs_infc_QPs)
			{
				constexpr auto ref_solution = 1.0;
				constexpr auto epsilon = 9.0E-12;
				EXPECT_NEAR(ref_solution, solution[0], epsilon);
				//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
			}
			for (const auto& solution : solution_at_ncs_infc_QPs)
			{
				constexpr auto ref_solution = 1.0;
				constexpr auto epsilon = 9.0E-12;
				EXPECT_NEAR(ref_solution, solution[0], epsilon);
				//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
			}
		}
	}
}
TEST(Discrete_Solution_DG, calculate_solution_at_post_points_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/RQ3.msh";	

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	constexpr auto post_order = 3;

	for (int solution_degree = 0; solution_degree < 8; ++solution_degree)
	{
		const auto governing_equation = std::make_shared<Burgers_2D>();
		const auto initial_condition = std::make_unique<Constant1>();
		auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		auto post_variable_convertor = std::make_unique<Node_Base_Convertor>(grid, post_order, *discrete_solution_DG);		
		
		const auto num_cells = grid.num_cells();
		for (int cell_index = 0; cell_index < num_cells; ++cell_index)
		{
			const auto solutions = discrete_solution_DG->calculate_solution_at_post_points(cell_index);

			for (const auto& solution : solutions)
			{
				constexpr auto ref_solution = 1.0;
				constexpr auto epsilon = 9.0E-12;
				EXPECT_NEAR(ref_solution, solution[0], epsilon);
				//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
			}
		}
	}
}
TEST(Discrete_Solution_DG, calculate_solution_at_post_element_centers_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/RQ3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	constexpr auto post_order = 3;

	for (int solution_degree = 0; solution_degree < 8; ++solution_degree)
	{
		const auto governing_equation = std::make_shared<Burgers_2D>();
		const auto initial_condition = std::make_unique<Constant1>();
		auto discrete_solution_DG = std::make_unique<Discrete_Solution_DG>(governing_equation, grid, *initial_condition, solution_degree);

		auto post_variable_convertor = std::make_unique<Center_Base_Convertor>(grid, post_order, *discrete_solution_DG);

		const auto num_cells = grid.num_cells();
		for (int cell_index = 0; cell_index < num_cells; ++cell_index)
		{
			const auto solutions = discrete_solution_DG->calculate_solution_at_post_element_centers(cell_index);

			for (const auto& solution : solutions)
			{
				constexpr auto ref_solution = 1.0;
				constexpr auto epsilon = 9.0E-12;
				EXPECT_NEAR(ref_solution, solution[0], epsilon);
				//EXPECT_DOUBLE_EQ(ref_solution, solution[0]); //suffer by numerical error
			}
		}
	}
}
TEST(Discrete_Solution_DG, turn_on_log_print)
{
	LOG << Log::print_on;
}