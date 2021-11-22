#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Governing_Equation.h"

#include <random>

TEST(Linear_Advection_2D, calculate_coordinate_projected_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_2D gov_eq({ 1.0,0.5 });

	std::vector<Euclidean_Vector> solutions(num);
	for (size_t i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}

	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num, { 1.0,0.5 });
	EXPECT_EQ(result, ref);
}
TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_2D gov_eq({ 1.0,0.5 });

	Euclidean_Vector solution_o = 1;
	Euclidean_Vector solution_n = 0;

	for (int i = 0; i < num; ++i) 
	{
		Euclidean_Vector normal = { 1.0 * i, 1.0 * i };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}
TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_2)
{
	constexpr size_t num = 5;

	Linear_Advection_2D gov_eq({ 1.0,0.5 });

	Euclidean_Vector solution_o = 1;
	Euclidean_Vector solution_n = 0;

	for (int i = 0; i < num; ++i) 
	{
		Euclidean_Vector normal = { -1.0 * i, -1.0 * i };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 1.5 * i;
		EXPECT_EQ(result, ref);
	}
}
TEST(Linear_Advection_2D, calculate_physical_flux_1)
{
	Linear_Advection_2D gov_eq({ 1.0,0.5 });

	std::vector<Euclidean_Vector> solutions = { {1},{2},{3},{4},{5} };

	for (int i = 0; i < solutions.size(); ++i)
	{
		const auto result = gov_eq.calculate_physical_flux(solutions[i]);
		const Matrix ref(1, 2, { 1.0 * (i + 1), 0.5 * (i + 1) });
		EXPECT_EQ(result, ref);
	}
}

TEST(Linear_Advection_3D, calculate_coordinate_projected_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_3D gov_eq({ 1.0, 0.5, 3.0 });

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}

	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num, { 1.0,0.5, 3.0 });
	EXPECT_EQ(result, ref);
}
TEST(Linear_Advection_3D, calculate_inner_face_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_3D gov_eq({ 1.0, 0.5, 2.0 });

	Euclidean_Vector solution_o = 1;
	Euclidean_Vector solution_n = 0;

	for (int i = 0; i < num; ++i) 
	{
		Euclidean_Vector normal = { 1.0 * i, 1.0 * i, 1.0 * i };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 3.5 * i;
		EXPECT_EQ(result, ref);
	}
}
TEST(Linear_Advection_3D, calculate_inner_face_maximum_lambdas_2)
{
	constexpr size_t num = 5;

	Linear_Advection_3D gov_eq({ 1.0, 0.5, 2.0 });

	Euclidean_Vector solution_o = 1;
	Euclidean_Vector solution_n = 0;

	for (int i = 0; i < num; ++i)
	{
		Euclidean_Vector normal = { -1.0 * i, -1.0 * i, -1.0 * i };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = 3.5 * i;
		EXPECT_EQ(result, ref);
	}
}
TEST(Linear_Advection_3D, calculate_physical_flux_1)
{
	Linear_Advection_3D gov_eq({ 1.0, 0.5, 2.0 });

	constexpr size_t num = 5;

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}

	for (int i = 0; i < solutions.size(); ++i)
	{
		const auto result = gov_eq.calculate_physical_flux(solutions[i]);
		const Matrix ref(1, 3, { 1.0 * (i), 0.5 * (i), 2.0 * i });
		EXPECT_EQ(result, ref);
	}
}

TEST(Burgers_2D, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr size_t num = 5;
	const Burgers_2D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (size_t i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}
	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num);
	for (int i = 0; i < num; ++i)
	{
		ref[i] = { 1.0 * i ,1.0 * i };
	}
	EXPECT_EQ(result, ref);
}
TEST(Burgers_2D, calculate_coordinate_projected_maximum_lambdas_2) {
	constexpr size_t num = 5;
	const Burgers_2D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
	{
		solutions[i] = { -1.0 * i };
	}
	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num);
	for (int i = 0; i < num; ++i)
	{
		ref[i] = { 1.0 * i ,1.0 * i };
	}
	EXPECT_EQ(result, ref);
}
TEST(Burgers_2D, calculate_inner_face_maximum_lambdas_1)
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	std::uniform_real_distribution<double> dis2(-1, 1);

	Burgers_2D gov_eq;

	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i)
	{
		Euclidean_Vector solution_o = { dis1(gen) };
		Euclidean_Vector solution_n = { dis1(gen) };

		Euclidean_Vector normal = { dis2(gen), dis2(gen) };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = std::max(std::abs(solution_o[0] * (normal[0] + normal[1])), std::abs(solution_n[0] * (normal[0] + normal[1])));
		EXPECT_EQ(result, ref);
	}
}
TEST(Burgers_2D, calculate_physical_flux_1) {
	constexpr size_t num = 5;

	const Burgers_2D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}

	for (int i = 0; i < num; ++i)
	{
		const auto result = gov_eq.calculate_physical_flux(solutions[i]);
		const Matrix ref = { 1,2,{ 0.5 * i * i, 0.5 * i * i } };
		EXPECT_EQ(result, ref);
	}
}

TEST(Burgers_3D, calculate_coordinate_projected_maximum_lambdas_1) {
	constexpr size_t num = 5;
	const Burgers_3D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (size_t i = 0; i < num; ++i)
		solutions[i] = { 1.0 * i };
	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num);
	for (int i = 0; i < num; ++i)
	{
		ref[i] = { 1.0 * i ,1.0 * i, 1.0 * i };
	}
	EXPECT_EQ(result, ref);
}
TEST(Burgers_3D, calculate_coordinate_projected_maximum_lambdas_2) 
{
	constexpr size_t num = 5;
	const Burgers_3D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
		solutions[i] = { -1.0 * i };
	const auto result = gov_eq.calculate_coordinate_projected_maximum_lambdas(solutions);

	std::vector<std::vector<double>> ref(num);
	for (int i = 0; i < num; ++i)
	{
		ref[i] = { 1.0 * i ,1.0 * i, 1.0 * i };
	}
	EXPECT_EQ(result, ref);
}
TEST(Burgers_3D, calculate_inner_face_maximum_lambdas_1) 
{
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	std::uniform_real_distribution<double> dis2(-1, 1);

	Burgers_3D gov_eq;
	constexpr size_t num = 5;
	for (size_t i = 0; i < num; ++i) 
	{
		Euclidean_Vector solution_o = { dis1(gen) };
		Euclidean_Vector solution_n = { dis1(gen) };

		Euclidean_Vector normal = { dis2(gen), dis2(gen), dis2(gen) };
		const auto result = gov_eq.calculate_inner_face_maximum_lambda(solution_o, solution_n, normal);

		const auto ref = std::max(std::abs(solution_o[0] * (normal[0] + normal[1] + normal[2])), std::abs(solution_n[0] * (normal[0] + normal[1] + normal[2])));
		EXPECT_EQ(result, ref);
	}
}
TEST(Burgers_3D, calculate_physical_flux_1) {
	constexpr size_t num = 5;

	const Burgers_3D gov_eq;

	std::vector<Euclidean_Vector> solutions(num);
	for (int i = 0; i < num; ++i)
	{
		solutions[i] = { 1.0 * i };
	}

	for (int i = 0; i < num; ++i)
	{
		const auto result = gov_eq.calculate_physical_flux(solutions[i]);
		const Matrix ref = { 1,3,{ 0.5 * i * i, 0.5 * i * i, 0.5 * i * i  } };
		EXPECT_EQ(result, ref);
	}
}




//TEST(Euler, conservative_to_primitive_1) {
//	constexpr ushort space_dimension = 2;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> solution = { 1,1,1,1 };
//	const auto result = Euler::conservative_to_primitive(solution);
//
//	const Euclidean_Vector<num_equation> ref = { 1,1,0,0 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Euler, conservative_to_primitive_2) {
//	constexpr ushort space_dimension = 3;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> solution = { 1,1,1,1,1.5 };
//	const auto result = Euler::conservative_to_primitive(solution);
//
//	const Euclidean_Vector<num_equation> ref = { 1,1,1,0,0 };
//	EXPECT_EQ(result, ref);
//}
//
//TEST(Euler, physical_flux_1) {
//	constexpr ushort space_dimension = 2;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> solution = { 1,1,1,1 };
//	const auto result = Euler::physical_flux(solution);
//
//	const Static_Matrix<num_equation,space_dimension> ref = { 1,1,1,1,1,1,1,1 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Euler, physical_flux_2) {
//	constexpr ushort space_dimension = 3;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> solution = { 1,1,1,1,1.5 };
//	const auto result = Euler::physical_flux(solution);
//
//	const Static_Matrix<num_equation, space_dimension> ref = { 1,1,1,1,1,1,1,1,1,1,1,1,1.5,1.5,1.5 };
//	EXPECT_EQ(result, ref);
//}
//
//TEST(Euler, coordinate_projected_maximum_lambda_1) {
//	constexpr ushort space_dimension = 2;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	std::vector<Euclidean_Vector<num_equation>> solution = { { 1,1,1,1 } };
//	const auto result = Euler::calculate_coordinate_projected_maximum_lambdas(solution);
//
//	const std::array<double, space_dimension> ref = { 1,1 };
//	EXPECT_EQ(result[0], ref);
//}
//TEST(Euler, coordinate_projected_maximum_lambda_2) {
//	constexpr ushort space_dimension = 3;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	std::vector<Euclidean_Vector<num_equation>> solution = { { 1,1,1,1,1.5 } };
//	const auto result = Euler::calculate_coordinate_projected_maximum_lambdas(solution);
//
//	const std::array<double, space_dimension> ref = { 1,1,1 };
//	EXPECT_EQ(result[0], ref);
//}
//
//TEST(Euler, inner_face_maximum_lambda_1) {
//	constexpr ushort space_dimension = 2;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> oc_solution = { 1,1,1,1 };
//	Euclidean_Vector<num_equation> nc_solution = { 1,0,0,0 };
//
//	Euclidean_Vector normal = { 1,0 };
//
//	const auto oc_pvariable = Euler::conservative_to_primitive(oc_solution);
//	const auto nc_pvariable = Euler::conservative_to_primitive(nc_solution);
//
//	const auto result = Euler::inner_face_maximum_lambda(oc_pvariable, nc_pvariable, normal);
//
//	const auto ref = 1.0;
//	EXPECT_EQ(result, ref);
//}
//TEST(Euler, inner_face_maximum_lambda_2) {
//	constexpr ushort space_dimension = 3;
//	constexpr ushort num_equation = space_dimension + 2;
//
//	Euclidean_Vector<num_equation> oc_solution = { 1,1,1,1,1.5 };
//	Euclidean_Vector<num_equation> nc_solution = { 1,0,0,0,0 };
//
//	Euclidean_Vector normal = { 1,0,0 };
//
//	const auto oc_pvariable = Euler::conservative_to_primitive(oc_solution);
//	const auto nc_pvariable = Euler::conservative_to_primitive(nc_solution);
//
//	const auto result = Euler::inner_face_maximum_lambda(oc_pvariable, nc_pvariable, normal);
//
//	const auto ref = 1.0;
//	EXPECT_EQ(result, ref);
//}