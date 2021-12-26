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

	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

	std::vector<std::vector<double>> ref(num, { 1.0,0.5 });
	EXPECT_EQ(result, ref);
}
TEST(Linear_Advection_2D, calculate_inner_face_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_2D gov_eq({ 1.0,0.5 });

	Euclidean_Vector solution_o = { 1 };
	Euclidean_Vector solution_n = { 0 };

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

	Euclidean_Vector solution_o = { 1 };
	Euclidean_Vector solution_n = { 0 };

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

	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

	std::vector<std::vector<double>> ref(num, { 1.0,0.5, 3.0 });
	EXPECT_EQ(result, ref);
}
TEST(Linear_Advection_3D, calculate_inner_face_maximum_lambdas_1)
{
	constexpr size_t num = 5;

	Linear_Advection_3D gov_eq({ 1.0, 0.5, 2.0 });

	Euclidean_Vector solution_o = { 1 };
	Euclidean_Vector solution_n = { 0 };

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

	Euclidean_Vector solution_o = { 1 };
	Euclidean_Vector solution_n = { 0 };

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
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

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
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

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
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

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
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solutions);

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


TEST(Euler_2D, extend_to_solution_1)
{
	Euler_2D gov_eq;
	Euclidean_Vector gov_sol = { 1,1,1,1 };
	gov_eq.extend_to_solution(gov_sol);

	const Euclidean_Vector ref = { 1,1,1,1,1,1,0,0 };
	EXPECT_EQ(gov_sol, ref);
}
TEST(Euler_2D, calculate_coordinate_projected_maximum_lambda_1)
{
	Euler_2D gov_eq;

	std::vector<Euclidean_Vector> solution = { { 1,1,1,1 } };
	gov_eq.extend_to_solution(solution[0]);
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solution);

	const std::vector<double> ref = { 1,1 };
	EXPECT_EQ(result[0], ref);
}
TEST(Euler_2D, calculate_inner_face_maximum_lambda_1) {
	Euler_2D gov_eq;

	Euclidean_Vector oc_solution = { 1,1,1,1 };
	Euclidean_Vector nc_solution = { 1,0,0,0 };
	Euclidean_Vector normal = { 1,0 };

	gov_eq.extend_to_solution(oc_solution);
	gov_eq.extend_to_solution(nc_solution);
	const auto result = gov_eq.calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal);

	const auto ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Euler_2D, calculate_physical_flux_1)
{
	Euler_2D gov_eq;

	Euclidean_Vector solution = { 1,1,1,1 };
	gov_eq.extend_to_solution(solution);
	const auto result = gov_eq.calculate_physical_flux(solution);

	const Matrix ref = { gov_eq.num_equations(), gov_eq.space_dimension(), { 1,1,1,1,1,1,1,1 } };
	EXPECT_EQ(result, ref);
}

TEST(Euler_3D, extend_to_solution_1)
{
	Euler_3D gov_eq;

	Euclidean_Vector gov_sol = { 1,1,1,1,1.5 };
	gov_eq.extend_to_solution(gov_sol);

	const Euclidean_Vector ref = { 1,1,1,1,1.5,1,1,1,0,0 };
	EXPECT_EQ(gov_sol, ref);
}
TEST(Euler_3D, calculate_coordinate_projected_maximum_lambda_1)
{
	Euler_3D gov_eq;

	std::vector<Euclidean_Vector> solution = { { 1,1,1,1,1.5 } };
	gov_eq.extend_to_solution(solution[0]);
	const auto result = gov_eq.calculate_cell_index_to_coordinate_projected_maximum_lambdas_table(solution);

	const std::vector<double> ref = { 1,1,1 };
	EXPECT_EQ(result[0], ref);
}
TEST(Euler_3D, calculate_inner_face_maximum_lambda_1)
{
	Euler_3D gov_eq;

	Euclidean_Vector oc_solution = { 1,1,1,1,1.5 };
	Euclidean_Vector nc_solution = { 1,0,0,0,0 };
	Euclidean_Vector normal = { 1,0,0 };

	gov_eq.extend_to_solution(oc_solution);
	gov_eq.extend_to_solution(nc_solution);
	const auto result = gov_eq.calculate_inner_face_maximum_lambda(oc_solution, nc_solution, normal);

	const auto ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Euler_3D, calculate_physical_flux_1)
{
	Euler_3D gov_eq;

	Euclidean_Vector solution = { 1,1,1,1,1.5 };
	gov_eq.extend_to_solution(solution);
	const auto result = gov_eq.calculate_physical_flux(solution);

	const Matrix ref = { gov_eq.num_equations(), gov_eq.space_dimension(), { 1,1,1,1,1,1,1,1,1,1,1,1,1.5,1.5,1.5 } };
	EXPECT_EQ(result, ref);
}