#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Governing_Equation.h"

#include <random>

GTEST_TEST(Linear_Advection_2D, calculate_physical_fluxes_1) {
	Linear_Advection_2D::initialize({ 1.0,0.5 });
	
	constexpr size_t num = 5;
	
	std::vector<Linear_Advection_2D::Solution_> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };	
	const auto result = Linear_Advection_2D::physical_fluxes(solutions);

	std::vector<Matrix<1, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { i, 0.5 * i };

	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection_2D, calculate_coordinate_projected_maximum_lambdas_1) {
	Linear_Advection_2D::initialize({ 1.0,0.5 });


	constexpr size_t num = 5;

	std::vector<Linear_Advection_2D::Solution_> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection_2D::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, 2>> ref(num, { 1.0,0.5 });
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_1) {	
	Linear_Advection_2D::initialize({ 1.0,0.5 });

	Linear_Advection_2D::Solution_ solution_o = 1;
	Linear_Advection_2D::Solution_ solution_n = 0;

	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i) {
		Linear_Advection_2D::Space_Vector_ normal = { i, i };
		const auto result = Linear_Advection_2D::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}
GTEST_TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_2) {
	Linear_Advection_2D::Solution_ solution_o;
	Linear_Advection_2D::Solution_ solution_n;

	constexpr size_t num = 5;
	for (int i = 0; i < num; ++i) {
		Linear_Advection_2D::Space_Vector_ normal = { -1 * i, -1 * i };
		const auto result = Linear_Advection_2D::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}

GTEST_TEST(Burgers_2D, calculate_physical_fluxes_1) {
	constexpr size_t num = 5;

	std::vector<Burgers_2D::Solution_> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers_2D::physical_fluxes(solutions);

	std::vector<Matrix<1, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { 0.5 * i * i, 0.5 * i * i };

	EXPECT_EQ(result, ref);
}

GTEST_TEST(Burgers_2D, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr size_t num = 5;

	std::vector<Burgers_2D::Solution_> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers_2D::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, 2>> ref(num);
	for (size_t i = 0; i < num; ++i) 
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Burgers_2D, calculate_coordinate_projected_maximum_lambdas_2) {
	constexpr size_t num = 5;

	std::vector<Burgers_2D::Solution_> solutions(num);
	for (int i = 0; i < num; ++i)
		solutions[i] = { -i };
	const auto result = Burgers_2D::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Burgers_2D, calculate_inner_face_maximum_lambdas_1) {	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	std::uniform_real_distribution<double> dis2(-1, 1);

	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i) {
		Burgers_2D::Solution_ solution_o = dis1(gen);
		Burgers_2D::Solution_ solution_n = dis1(gen);

		Burgers_2D::Space_Vector_ normal = { dis2(gen), dis2(gen) };
		const auto result = Burgers_2D::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = std::max(std::abs(solution_o[0] * (normal[0]+normal[1])), std::abs(solution_n[0] * (normal[0] + normal[1])));
		EXPECT_EQ(result, ref);
	}
}