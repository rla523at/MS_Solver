#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Reconstruction_Method.h"

TEST(Reconstruction_Method, activation_function_1) {
	HardTanh f;

	Dynamic_Euclidean_Vector v = { -1,0,0.5,1,2,3 };
	v.apply(f);

	Dynamic_Euclidean_Vector ref = { 0,0,0.5,1,1,1 };
	EXPECT_EQ(v, ref);
}
TEST(Reconstruction_Method, activation_function_2) {
	ReLU f;

	Dynamic_Euclidean_Vector v = { -1,0,0.5,1,2,3 };
	v.apply(f);

	Dynamic_Euclidean_Vector ref = { 0,0,0.5,1,2,3 };
	EXPECT_EQ(v, ref);
}


TEST(hMLP_Reconstruction, is_satisfy_P1_projected_MLP_condition_1) {

	constexpr ushort num_equation = 1;
	constexpr ushort space_dimension = 2;
	constexpr ushort polynomial_order = 1;

	const auto P1_projected_solution = 1.5;
	const auto allowable_min = 0.5;
	const auto allowable_max = 2.0;

	const auto result = hMLP_Reconstruction<num_equation, space_dimension, polynomial_order>::is_satisfy_P1_projected_MLP_condition(P1_projected_solution, allowable_min, allowable_max);
	const auto ref = true;
	EXPECT_EQ(result, ref);
}
TEST(hMLP_Reconstruction, is_satisfy_P1_projected_MLP_condition_2) {

	constexpr ushort num_equation = 1;
	constexpr ushort space_dimension = 2;
	constexpr ushort polynomial_order = 1;

	const auto P1_projected_solution = 2.5;
	const auto allowable_min = 0.5;
	const auto allowable_max = 2.0;

	const auto result = hMLP_Reconstruction<num_equation, space_dimension, polynomial_order>::is_satisfy_P1_projected_MLP_condition(P1_projected_solution, allowable_min, allowable_max);
	const auto ref = false;
	EXPECT_EQ(result, ref);
}