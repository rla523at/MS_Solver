#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Euclidean_Vector.h"

GTEST_TEST(Euclidean_Vector, constructor_1) {
	const Euclidean_Vector<2> result = { 1 };

	const Euclidean_Vector ref = { 1,0 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Euclidean_Vector, operator_addition_1) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 5,7,9 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Euclidean_Vector, operator_addition_2) {
	const Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = { 4,5,6,7 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 5,7,9,11 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Euclidean_Vector, operator_substraction_1) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 - v2;

	const Euclidean_Vector ref = { -3,-3,-3 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Euclidean_Vector, operator_scalar_multiplication_1) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const auto result = v1 * 2;

	const Euclidean_Vector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Euclidean_Vector, operator_scalar_multiplication_2) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const auto result = 2 * v1;

	const Euclidean_Vector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Euclidean_Vector, inner_product_1) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1.inner_product(v2);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Euclidean_Vector, norm_1) {
	const Euclidean_Vector v1 = { 3,4 };
	const auto result = v1.L2_norm();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}

TEST(Dynamic_Euclidean_Vector, be_absolute_1) {
	Dynamic_Euclidean_Vector v1 = { -1, -2 ,3 ,-4 ,5 };
	v1.be_absolute();

	const Dynamic_Euclidean_Vector ref = { 1,2,3,4,5 };
	EXPECT_EQ(v1, ref);
}

TEST(Dynamic_Euclidean_Vector, inner_product_1) {
	Dynamic_Euclidean_Vector v1 = { 3,4 };
	const auto result = v1.inner_product(v1);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}

TEST(ms, min_value_gathering_vector_1) {
	Euclidean_Vector v1 = { 1,2,3 };
	Euclidean_Vector v2 = { 2,3,1 };
	Euclidean_Vector v3 = { 3,1,2 };

	std::vector<Euclidean_Vector<3>> vec = { v1,v2,v3 };

	const auto result = ms::gather_min_value(vec);
	const Euclidean_Vector ref = { 1,1,1 };
	EXPECT_EQ(result, ref);
}

TEST(ms, max_value_gathering_vector_1) {
	Euclidean_Vector v1 = { 1,2,3 };
	Euclidean_Vector v2 = { 2,3,1 };
	Euclidean_Vector v3 = { 3,2,1 };

	std::vector<Euclidean_Vector<3>> vec = { v1,v2,v3 };

	const auto result = ms::gather_max_value(vec);
	const Euclidean_Vector ref = { 3,3,3 };
	EXPECT_EQ(result, ref);
}