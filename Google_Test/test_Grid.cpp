#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Grid.h"

TEST(Grid, turn_off_log_print)
{
	LOG << Log::print_off;
}
TEST(Grid, num_pbdry_pairs_1) 
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/Quad3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);
	
	const auto result = grid.num_periodic_boundary_pairs();

	const auto ref = 6;
	EXPECT_EQ(ref, result);
}
TEST(Grid, pbdry_oc_nc_index_pair_1)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/Quad3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	const auto result = grid.periodic_boundary_oc_nc_index_pair(0);

	const std::pair<uint, uint> ref = { 2,0 };
	EXPECT_EQ(ref, result);
}
TEST(Grid, pbdry_oc_nc_index_pair_2)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/Quad3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	const auto result = grid.periodic_boundary_oc_nc_index_pair(3);

	const std::pair<uint, uint> ref = { 0,6 };
	EXPECT_EQ(ref, result);
}
TEST(Grid, pbdry_oc_nc_index_pair_3)
{
	constexpr auto space_dimension = 2;
	constexpr auto grid_file_path = "RSC/Grid/2D/Quad3.msh";

	Gmsh_Convertor convertor(space_dimension);
	Grid grid(convertor, grid_file_path);

	const auto result = grid.periodic_boundary_oc_nc_index_pair(5);

	const std::pair<uint, uint> ref = { 2,8 };
	EXPECT_EQ(ref, result);
}
TEST(Grid, turn_on_log_print)
{
	LOG << Log::print_on;
}

//
//TEST(Grid, vnode_index_to_matched_vnode_index_set_1) 
//{
//	constexpr ushort space_dimension = 2;
//
//	auto grid_elements = Grid_Element_Builder<Gmsh, space_dimension>::build_from_grid_file("Quad3");
//	const auto grid = Grid<space_dimension>(std::move(grid_elements));
//	const auto result = grid.pbdry_vnode_index_to_matched_pbdry_vnode_index_set().at(0);
//
//	const std::set<uint> ref = { 3,12,15 };
//	EXPECT_EQ(ref, result);
//}
//
//TEST(Grid, vnode_index_to_share_cell_indexes_1) {
//	constexpr ushort space_dimension = 2;
//
//	auto grid_elements = Grid_Element_Builder<Gmsh, space_dimension>::build_from_grid_file("Quad3");
//	const auto grid = Grid<space_dimension>(std::move(grid_elements));
//	const auto result = grid.get_vnode_index_to_share_cell_index_set_consider_pbdry().at(0);
//
//	const std::set<uint> ref = { 0,2,6,8 };
//	EXPECT_EQ(ref, result);
//}

//
//GTEST_TEST(Grid_Info_Extractor, volume) {
//	auto grid_data = Grid_File_Convertor<Gmsh, 2>::convert("RSC/Grid/Quad_10.msh");
//	const auto grid_info = Grid_Info_Extractor<2>::convert(std::move(grid_data));
//
//	const auto& cell_volumes = grid_info.cell_grid_information.volumes;
//	
//	const auto ref = 0.01;
//	for (const auto& volume : cell_volumes)
//		EXPECT_NEAR(volume, ref, 1.0E-16);
//		//EXPECT_DOUBLE_EQ(volume, ref);	//suffer by round off error
//}
//
//GTEST_TEST(Grid_Info_Extractor, coordinate_projected_volumes) {
//	auto grid_data = Grid_File_Convertor<Gmsh, 2>::convert("RSC/Grid/Quad_10.msh");
//	const auto grid_info = Grid_Info_Extractor<2>::convert(std::move(grid_data));
//
//	const auto& coordinate_projected_volumes = grid_info.cell_grid_information.coordinate_projected_volumes;
//
//	const auto ref = 0.1;
//	for (const auto& coordinate_projected_volumes : coordinate_projected_volumes) {
//		const auto [x_projected, y_projected] = coordinate_projected_volumes;
//		EXPECT_NEAR(x_projected, ref, 9.0E-16);
//		EXPECT_NEAR(y_projected, ref, 9.0E-16);	//suffer by round off error
//	}
//}