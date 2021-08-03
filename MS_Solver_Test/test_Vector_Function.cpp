#include "gtest/gtest.h"
#include "../MS_Solver/INC/Vector_Function.h"
#include "../MS_Solver/INC/Polynomial.h"

GTEST_TEST(Vector_Function, range_dimension)
{
	Vector_Function<Polynomial<2>> vf(2);
	const auto result = vf.range_dimension();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Vector_Function, operator_access_1)
{
	Polynomial<2> x("x0");
	Polynomial<2> y("x1");

	Vector_Function<Polynomial<2>> vf(2);
	vf[1] = 2 * x + y;

	const auto ref = 2 * x + y;
	EXPECT_EQ(ref, vf[1]);
}

GTEST_TEST(Vector_Function, operator_calculate_1)
{
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>> vf = { x + y , 2 * x + y };	
	Euclidean_Vector v = { 1,1 };
	const auto result = vf(v);
	
	Dynamic_Euclidean_Vector_ ref = { 2,3 };
	EXPECT_EQ(ref, result);
}

GTEST_TEST(Vector_Function, differentiate_1)
{
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>> vf = { x + y , 2 * x + y };
	const auto result = vf.differentiate<0>();

	Vector_Function<Polynomial<domain_dimension>> ref = { 1,2 };
	EXPECT_EQ(ref, result);
}