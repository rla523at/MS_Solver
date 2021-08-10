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

TEST(Matrix, column_1) {
	const Matrix<2, 3> m = { 1,2,3,4,5,6 };
	const auto result = m.column(0);

	const Euclidean_Vector ref = { 1,4 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, column_2) {
	const Matrix<2, 3> m = { 1,2,3,4,5,6 };
	const auto result = m.column(2);

	const Euclidean_Vector ref = { 3,6 };
	EXPECT_EQ(result, ref);
}

TEST(Dynamic_Matrix, change_rows_1) {
	const Matrix<2, 2> m = { 1,2,3,4 };

	Dynamic_Matrix_ dm(4, 2);
	dm.change_rows(0, m);

	const Dynamic_Matrix_ ref(4, 2, { 1,2,3,4,0,0,0,0 });
	EXPECT_EQ(dm, ref);
}

TEST(Dynamic_Matrix, change_columns_1) {
	const Matrix<2, 2> m = { 1,2,3,4 };

	Dynamic_Matrix_ dm(2, 4);
	dm.change_columns(0, m);

	const Dynamic_Matrix_ ref(2, 4, { 1,2,0,0,3,4,0,0 });
	EXPECT_EQ(dm, ref);
}