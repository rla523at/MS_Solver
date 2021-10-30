#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Governing_Equation.h"

#include <random>

GTEST_TEST(Linear_Advection, calculate_physical_fluxes_1) {
	constexpr ushort space_dimension = 2;

	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	constexpr size_t num = 5;
	
	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };	
	const auto result = Linear_Advection<space_dimension>::physical_fluxes(solutions);

	std::vector<Static_Matrix<1, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { i, 0.5 * i };

	EXPECT_EQ(result, ref);
}
GTEST_TEST(Linear_Advection, calculate_physical_fluxes_2) {
	constexpr ushort space_dimension = 3;

	Linear_Advection<space_dimension>::initialize({ 1.0, 0.5, 2.0 });

	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection<space_dimension>::physical_fluxes(solutions);

	std::vector<Static_Matrix<1, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { i, 0.5 * i, 2.0 * i };

	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, 2>> ref(num, { 1.0,0.5 });
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Linear_Advection, calculate_coordinate_projected_maximum_lambdas_2) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0, 0.5, 3.0 });

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Linear_Advection<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, space_dimension>> ref(num, { 1.0,0.5, 3.0 });
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Linear_Advection, calculate_inner_face_maximum_lambdas_1) {	
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	Euclidean_Vector<1> solution_o = 1;
	Euclidean_Vector<1> solution_n = 0;

	for (size_t i = 0; i < num; ++i) {
		Euclidean_Vector<space_dimension> normal = { i, i };
		const auto result = Linear_Advection<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}
GTEST_TEST(Linear_Advection, calculate_inner_face_maximum_lambdas_2) {
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	Euclidean_Vector<1> solution_o = 1;
	Euclidean_Vector<1> solution_n = 0;

	for (int i = 0; i < num; ++i) {
		Euclidean_Vector<space_dimension> normal = { -i, -i };
		const auto result = Linear_Advection<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}
GTEST_TEST(Linear_Advection, calculate_inner_face_maximum_lambdas_3) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0, 0.5, 2.0 });

	Euclidean_Vector<1> solution_o = 1;
	Euclidean_Vector<1> solution_n = 0;

	
	for (size_t i = 0; i < num; ++i) {
		Euclidean_Vector<space_dimension> normal = { i, i, i };
		const auto result = Linear_Advection<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 3.5 * i;
		EXPECT_EQ(result, ref);
	}
}
GTEST_TEST(Linear_Advection, calculate_inner_face_maximum_lambdas_4) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	Linear_Advection<space_dimension>::initialize({ 1.0, 0.5, 2.0 });

	Euclidean_Vector<1> solution_o = 1;
	Euclidean_Vector<1> solution_n = 0;

	for (int i = 0; i < num; ++i) {
		Euclidean_Vector<space_dimension> normal = { -i, -i, -i };
		const auto result = Linear_Advection<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 3.5 * i;
		EXPECT_EQ(result, ref);
	}
}


GTEST_TEST(Burgers, calculate_physical_fluxes_1) {
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers<space_dimension>::physical_fluxes(solutions);

	std::vector<Static_Matrix<1, 2>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { 0.5 * i * i, 0.5 * i * i };

	EXPECT_EQ(result, ref);
}
GTEST_TEST(Burgers, calculate_physical_fluxes_2) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers<space_dimension>::physical_fluxes(solutions);

	std::vector<Static_Matrix<1, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { 0.5 * i * i, 0.5 * i * i, 0.5 * i * i };

	EXPECT_EQ(result, ref);
}

GTEST_TEST(Burgers, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i) 
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Burgers, calculate_coordinate_projected_maximum_lambdas_2) {
	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (int i = 0; i < num; ++i)
		solutions[i] = { -i };
	const auto result = Burgers<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Burgers, calculate_coordinate_projected_maximum_lambdas_3) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { i };
	const auto result = Burgers<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i) ,static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Burgers, calculate_coordinate_projected_maximum_lambdas_4) {
	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;

	std::vector<Euclidean_Vector<1>> solutions(num);
	for (int i = 0; i < num; ++i)
		solutions[i] = { -i };
	const auto result = Burgers<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::array<double, space_dimension>> ref(num);
	for (size_t i = 0; i < num; ++i)
		ref[i] = { static_cast<double>(i) ,static_cast<double>(i),static_cast<double>(i) };
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Burgers, calculate_inner_face_maximum_lambdas_1) {	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	std::uniform_real_distribution<double> dis2(-1, 1);

	constexpr ushort space_dimension = 2;
	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i) {
		Euclidean_Vector<1> solution_o = dis1(gen);
		Euclidean_Vector<1> solution_n = dis1(gen);

		Euclidean_Vector<space_dimension> normal = { dis2(gen), dis2(gen) };
		const auto result = Burgers<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = std::max(std::abs(solution_o[0] * (normal[0]+normal[1])), std::abs(solution_n[0] * (normal[0] + normal[1])));
		EXPECT_EQ(result, ref);
	}
}
GTEST_TEST(Burgers, calculate_inner_face_maximum_lambdas_2) {
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	std::uniform_real_distribution<double> dis2(-1, 1);

	constexpr ushort space_dimension = 3;
	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i) {
		Euclidean_Vector<1> solution_o = dis1(gen);
		Euclidean_Vector<1> solution_n = dis1(gen);

		Euclidean_Vector<space_dimension> normal = { dis2(gen), dis2(gen), dis2(gen) };
		const auto result = Burgers<space_dimension>::inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = std::max(std::abs(solution_o[0] * (normal[0] + normal[1]+ normal[2])), std::abs(solution_n[0] * (normal[0] + normal[1] + normal[2])));
		EXPECT_EQ(result, ref);
	}
}

TEST(Euler, conservative_to_primitive_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> solution = { 1,1,1,1 };
	const auto result = Euler<space_dimension>::conservative_to_primitive(solution);

	const Euclidean_Vector<num_equation> ref = { 1,1,0,0 };
	EXPECT_EQ(result, ref);
}
TEST(Euler, conservative_to_primitive_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> solution = { 1,1,1,1,1.5 };
	const auto result = Euler<space_dimension>::conservative_to_primitive(solution);

	const Euclidean_Vector<num_equation> ref = { 1,1,1,0,0 };
	EXPECT_EQ(result, ref);
}

TEST(Euler, physical_flux_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> solution = { 1,1,1,1 };
	const auto result = Euler<space_dimension>::physical_flux(solution);

	const Static_Matrix<num_equation,space_dimension> ref = { 1,1,1,1,1,1,1,1 };
	EXPECT_EQ(result, ref);
}
TEST(Euler, physical_flux_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> solution = { 1,1,1,1,1.5 };
	const auto result = Euler<space_dimension>::physical_flux(solution);

	const Static_Matrix<num_equation, space_dimension> ref = { 1,1,1,1,1,1,1,1,1,1,1,1,1.5,1.5,1.5 };
	EXPECT_EQ(result, ref);
}

TEST(Euler, coordinate_projected_maximum_lambda_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_equation = space_dimension + 2;

	std::vector<Euclidean_Vector<num_equation>> solution = { { 1,1,1,1 } };
	const auto result = Euler<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solution);

	const std::array<double, space_dimension> ref = { 1,1 };
	EXPECT_EQ(result[0], ref);
}
TEST(Euler, coordinate_projected_maximum_lambda_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_equation = space_dimension + 2;

	std::vector<Euclidean_Vector<num_equation>> solution = { { 1,1,1,1,1.5 } };
	const auto result = Euler<space_dimension>::calculate_coordinate_projected_maximum_lambdas(solution);

	const std::array<double, space_dimension> ref = { 1,1,1 };
	EXPECT_EQ(result[0], ref);
}

TEST(Euler, inner_face_maximum_lambda_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> oc_solution = { 1,1,1,1 };
	Euclidean_Vector<num_equation> nc_solution = { 1,0,0,0 };

	Euclidean_Vector<space_dimension> normal = { 1,0 };

	const auto oc_pvariable = Euler<space_dimension>::conservative_to_primitive(oc_solution);
	const auto nc_pvariable = Euler<space_dimension>::conservative_to_primitive(nc_solution);

	const auto result = Euler<space_dimension>::inner_face_maximum_lambda(oc_pvariable, nc_pvariable, normal);

	const auto ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Euler, inner_face_maximum_lambda_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_equation = space_dimension + 2;

	Euclidean_Vector<num_equation> oc_solution = { 1,1,1,1,1.5 };
	Euclidean_Vector<num_equation> nc_solution = { 1,0,0,0,0 };

	Euclidean_Vector<space_dimension> normal = { 1,0,0 };

	const auto oc_pvariable = Euler<space_dimension>::conservative_to_primitive(oc_solution);
	const auto nc_pvariable = Euler<space_dimension>::conservative_to_primitive(nc_solution);

	const auto result = Euler<space_dimension>::inner_face_maximum_lambda(oc_pvariable, nc_pvariable, normal);

	const auto ref = 1.0;
	EXPECT_EQ(result, ref);
}