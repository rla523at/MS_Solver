#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Matrix.h"

GTEST_TEST(Matrix, operator_mv_1) {
	const Euclidean_Vector v = { 1,2 };
	const Matrix<2,2> m = { 1,2,3,4 };
	const auto result = m * v;

	const Euclidean_Vector ref = { 5,11 };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Matrix, operator_mv_2) {
	const Matrix<1, 2> m = { 1,1 };
	const Euclidean_Vector v = { 1,2 };
	
	const auto result = (m+m) * v;

	const Euclidean_Vector ref = 6;
	EXPECT_EQ(result, ref);
}

TEST(Matrix, scalar_multiplication_1) {
	std::array<double, 9> ar;
	ar.fill(1);
	
	const Matrix<3, 3> m(ar);
	const auto result = m * 10;

	std::array<double, 9> ref_val;
	ref_val.fill(10);

	const Matrix<3, 3> ref(ref_val);
	EXPECT_EQ(result, ref);
}
TEST(Matrix, scalar_multiplication_2) {
	std::array<double, 900> ar;
	ar.fill(1);

	const Matrix<30, 30> m(ar);
	const auto result = m * 10;

	std::array<double, 900> ref_val;
	ref_val.fill(10);

	const Matrix<30, 30> ref(ref_val);
	EXPECT_EQ(result, ref);
}