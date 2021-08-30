#pragma once
#include "gtest/gtest.h"

#include "../MS_Solver/INC/Grid_Builder.h"


TEST(Grid, set_of_face_share_cell_indexes) {
	constexpr ushort space_dimension = 2;
	
	const auto grid = Grid_Builder<space_dimension>::build<Gmsh>("Quad3");
	const auto set_of_face_share_cell_indexes = grid.calculate_set_of_face_share_cell_indexes();

	auto result = set_of_face_share_cell_indexes[0];
	std::sort(result.begin(), result.end());

	const std::vector<size_t> ref = { 1,2,3,6 };
	EXPECT_EQ(ref, result);
}

TEST(Grid, vnode_index_to_matched_vnode_index_set_1) {
	constexpr ushort space_dimension = 2;

	const auto grid = Grid_Builder<space_dimension>::build<Gmsh>("Quad3");
	const auto result = grid.connectivity.vnode_index_to_matched_vnode_index_set.at(0);

	const std::set<uint> ref = { 3,12,15 };
	EXPECT_EQ(ref, result);
}

TEST(Grid, vnode_index_to_share_cell_indexes_1) {
	constexpr ushort space_dimension = 2;

	const auto grid = Grid_Builder<space_dimension>::build<Gmsh>("Quad3");
	const auto result = grid.connectivity.vnode_index_to_share_cell_index_set.at(0);

	const std::set<uint> ref = { 0,2,6,8 };
	EXPECT_EQ(ref, result);
}

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