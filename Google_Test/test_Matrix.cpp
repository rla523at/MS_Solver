#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Matrix.h"
#include "../MS_Solver/INC/Euclidean_Vector.h"

#include <array>

TEST(Matrix, construct_diagonal_matrix1) {
	Matrix m(2, {3, 3});
	Matrix ref(2, 2, { 3,0,0,3 });

	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_row1) {
	Matrix m(2, 2, { 1,2,3,4 });

	std::vector<double> v = { 5,6 };
	m.change_row(0, v);

	Matrix ref(2, 2, {5, 6, 3, 4});
	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_row2) {
	Matrix m(2, 2, { 1,2,3,4 });

	std::array<double,2> v = { 5,6 };
	m.change_row(0, v);

	Matrix ref(2, 2, { 5, 6, 3, 4 });
	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_rows1) {
	const Matrix m = { 2, 2,  {1,2,3,4} };

	Matrix dm(4, 2);
	dm.change_rows(0, m);

	const Matrix ref(4, 2, { 1,2,3,4,0,0,0,0 });
	EXPECT_EQ(dm, ref);
}
TEST(Matrix, operator_multiplication1) {
	const Matrix m1 = { 1,2, {1,2} };
	const Matrix m2 = { 2,2,{ 2,2,2,2 } };
	const auto result = m1 * m2;

	const Matrix ref = { 1,2,{ 6,6 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_multiplication2) {
	const Matrix m1 = { 3,{1,0,0} };
	const Matrix m2 = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const auto result = m1 * m2;

	const Matrix ref = { 3,3,{ 1,2,3,0,0,0,0,0,0 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_multiplication3) {
	const Matrix m1 = { 3,{1,0,0} };
	const Matrix m2 = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const auto result = m2 * m1;
	const Matrix ref = { 3,3,{ 1,0,0,4,0,0,7,0,0 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_multiplication4)
{
	const Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const Euclidean_Vector v = { -1,-1,1 };
	const auto result = m * v;
	const Euclidean_Vector ref = { 0,-3,-6 };
	EXPECT_EQ(result, ref);
}

TEST(Matrix, row1) {
	Matrix dm(2, 3, { 1,2,3,4,5,6 });
	const auto result = dm.row(0);

	const std::vector<double> ref = { 1,2,3 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, be_inverse_1) {
	Matrix m(2, 2, { 2,1,1,1 });
	m.be_inverse();

	Matrix ref(2, 2, { 1,-1,-1,2 });
	EXPECT_EQ(m, ref);
}
TEST(Matrix, inverse_1) {
	Matrix m(2, 2, { 2,1,1,1 });

	const auto result = m.inverse();	

	Matrix ref(2, 2, { 1,-1,-1,2 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, change_columns_1)
{
	const Matrix m = { 2,2,{ 1,2,3,4 } };

	Matrix dm(2, 4);
	dm.change_columns(0, m);

	const Matrix ref(2, 4, { 1,2,0,0,3,4,0,0 });
	EXPECT_EQ(dm, ref);
}
TEST(Matrix, change_columns_2)
{
	Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
	const Matrix m1 = { 2,2,{ 1,1,4,4 } };
	m.change_columns(1, m1);

	const Matrix ref = { 2,3,{ 1,1,1,4,4,4 } };
	EXPECT_EQ(m, ref);
}



//TEST(Matrix, construct_diagonal_matrix_1) {
//	std::array<double, 2> ar = { 2,2 };
//	Matrix<2, 2> result = ar;
//
//	Matrix<2, 2> ref = { 2,0,0,2 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, construct_diagonal_matrix_2) {
//	std::array<double, 2> ar = { 2,2 };
//	Matrix result = ar;
//
//	Matrix<2, 2> ref = { 2,0,0,2 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_plus_1) {
//	std::array<double, 16> ar;
//	ar.fill(1);
//
//	const Matrix<4, 4> m1 = ar;
//	const Matrix<4, 4> m2 = ar;
//	const auto result = m1 + m2;
//
//	std::array<double, 16> ref_val;
//	ref_val.fill(2);
//
//	const Matrix<4, 4> ref = ref_val;
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_plus_2) {
//	std::array<double, 36> ar;
//	ar.fill(1);
//
//	const Matrix<6, 6> m1 = ar;
//	const Matrix<6, 6> m2 = ar;
//	const auto result = m1 + m2;
//
//	std::array<double, 36> ref_val;
//	ref_val.fill(2);
//
//	const Matrix<6, 6> ref = ref_val;
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_mv_1) {
//	const Euclidean_Vector v = { 1,2 };
//	const Matrix<2,2> m = { 1,2,3,4 };
//	const auto result = m * v;
//
//	const Euclidean_Vector ref = { 5,11 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, operator_mv_2) {
//	const Matrix<1, 2> m = { 1,1 };
//	const Euclidean_Vector v = { 1,2 };
//	
//	const auto result = (m+m) * v;
//
//	const Euclidean_Vector ref = 6;
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, scalar_multiplication_1) {
//	std::array<double, 9> ar;
//	ar.fill(1);
//	
//	const Matrix<3, 3> m(ar);
//	const auto result = m * 10;
//
//	std::array<double, 9> ref_val;
//	ref_val.fill(10);
//
//	const Matrix<3, 3> ref(ref_val);
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, scalar_multiplication_2) {
//	std::array<double, 100> ar;
//	ar.fill(1);
//
//	const Matrix<10, 10> m(ar);
//	const auto result = m * 10;
//
//	std::array<double, 100> ref_val;
//	ref_val.fill(10);
//
//	const Matrix<10, 10> ref(ref_val);
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, gemm_1) {
//	Matrix<1, 2> m = { 1,2 };
//	const Matrix<2, 2> m2 = { 2,2,2,2 };
//	m *= m2;
//
//	const Matrix<1, 2> ref = { 6,6 };
//	EXPECT_EQ(m, ref);
//}
//TEST(Matrix, gemm_2) {
//	std::array<double, 3> ar = { 1, 0, 0 };
//	Matrix m = ar;
//
//	const Matrix<3, 3> m2 = { 1,2,3,4,5,6,7,8,9 };
//	const auto result = m * m2;
//
//	const Matrix<3, 3> ref = { 1,2,3,0,0,0,0,0,0 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, gemm_3) {
//	std::array<double, 3> ar = { 1, 0, 0 };
//	Matrix m = ar;
//
//	const Matrix<3, 3> m2 = { 1,2,3,4,5,6,7,8,9 };
//	const auto result = m2 * m;
//
//	const Matrix<3, 3> ref = { 1,0,0,4,0,0,7,0,0 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, column_1) {
//	const Matrix<2, 3> m = { 1,2,3,4,5,6 };
//	const auto result = m.column(0);
//
//	const Euclidean_Vector ref = { 1,4 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Matrix, column_2) {
//	const Matrix<2, 3> m = { 1,2,3,4,5,6 };
//	const auto result = m.column(2);
//
//	const Euclidean_Vector ref = { 3,6 };
//	EXPECT_EQ(result, ref);
//}

