#include "gtest/gtest.h"
#include "../MS_Solver/INC/Vector_Function.h"

GTEST_TEST(Vector_Function, range_dimension)
{
	Vector_Function<2> vf(2);
	const auto result = vf.range_dimension();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Vector_Function, operator_access_1)
{
	Polynomial<2> x("x0");
	Polynomial<2> y("x1");

	Vector_Function<2> vf(2);
	vf[1] = 2 * x + y;

	const auto ref = 2 * x + y;
	EXPECT_EQ(ref, vf[1]);
}