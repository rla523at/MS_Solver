#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Boundary_Flux_Function.h"

GTEST_TEST(Supersonic_Outlet, Linear_Advection_1) {
	constexpr ushort space_dimension = 2;
	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });
	
	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Linear_Advection<space_dimension>>>::make(boundary_type);				

	const Euclidean_Vector solution = { 2 };
	const Euclidean_Vector normal = { 1,0 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 2 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet, Linear_Advection_2) {
	constexpr ushort space_dimension = 2;
	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Linear_Advection<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector solution = { 2 };
	const Euclidean_Vector normal = { 1,1 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 3 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet, Linear_Advection_3) {
	constexpr ushort space_dimension = 2;
	Linear_Advection<space_dimension>::initialize({ 1.0,0.5 });

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Linear_Advection<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector solution = { 2 };
	const Euclidean_Vector normal = { 1,-1 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 1 };
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Supersonic_Outlet, Burgers_1) {
	constexpr ushort space_dimension = 2;

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Burgers<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector solution = { 4 };
	const Euclidean_Vector normal = { 1,0 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 8 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet, Burgers_2) {
	constexpr ushort space_dimension = 2;

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Burgers<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector solution = { 4 };
	const Euclidean_Vector normal = { 1,1 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 16 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Supersonic_Outlet, Burgers_3) {
	constexpr ushort space_dimension = 2;

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Burgers<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector solution = { 4 };
	const Euclidean_Vector normal = { 1,-1 };
	const auto result = bff->calculate(solution, normal);

	const Euclidean_Vector ref = { 0 };
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Supersonic_Outlet, Euler_1) {
	constexpr ushort space_dimension = 2;

	const auto boundary_type = ElementType::supersonic_outlet;
	const auto bff = Boundary_Flux_Function_Factory<LLF<Euler<space_dimension>>>::make(boundary_type);

	const Euclidean_Vector cvariable = { 1,0,0,1 };
	const Euclidean_Vector normal = { 1,0 };
	const auto result = bff->calculate(cvariable, normal);

	const Euclidean_Vector ref = { 0, (1.4 - 1), 0, 0 };
	EXPECT_EQ(result, ref);
}
