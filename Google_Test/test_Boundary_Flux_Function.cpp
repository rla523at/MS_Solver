#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Boundary_Flux_Function.h"

TEST(Initial_Constant_BC, Linear_Advection_2D_1) 
{
	const auto gov_eq = std::make_shared<Linear_Advection_2D>(1.0, 0.5);
	const auto numerical_flux_function = std::make_shared<LLF>(gov_eq);
	Initial_Constant_BC boundary_flux_function(numerical_flux_function);

	for (int i = 0; i < 10; ++i)
	{
		const Euclidean_Vector solution = { 1.0 * i + 1.0 };
		const Euclidean_Vector normal = { 1,0 };
		const auto result = boundary_flux_function.calculate(solution, normal);

		const Euclidean_Vector ref = { 1.0 + i };
		EXPECT_EQ(result, ref);
	}
}
TEST(Initial_Constant_BC, Linear_Advection_2D_2) 
{
	const auto gov_eq = std::make_shared<Linear_Advection_2D>(1.5, 0.5);
	const auto numerical_flux_function = std::make_shared<LLF>(gov_eq);
	Initial_Constant_BC boundary_flux_function(numerical_flux_function);

	for (int i = 0; i < 10; ++i)
	{
		const Euclidean_Vector solution = { 1.0 * i + 1.0 };
		const Euclidean_Vector normal = { std::sqrt(2) / 2.0, std::sqrt(2) / 2.0 };
		const auto result = boundary_flux_function.calculate(solution, normal);

		const Euclidean_Vector ref = { (1.0 + i) *std::sqrt(2) };
		EXPECT_DOUBLE_EQ(result[0], ref[0]);
	}
}

TEST(Initial_Constant_BC, Burgers_2D_1) 
{
	const auto gov_eq = std::make_shared<Burgers_2D>();
	const auto numerical_flux_function = std::make_shared<LLF>(gov_eq);
	Initial_Constant_BC boundary_flux_function(numerical_flux_function);

	for (int i = 0; i < 10; ++i)
	{
		const Euclidean_Vector solution = { 1.0 * i + 1.0 };
		const Euclidean_Vector normal = { 1,0 };
		const auto result = boundary_flux_function.calculate(solution, normal);

		const Euclidean_Vector ref = { 0.25 * (normal[0] + normal[1]) * (3 * i * i + 4 * i + 2) };
		EXPECT_EQ(result, ref);
	}
}