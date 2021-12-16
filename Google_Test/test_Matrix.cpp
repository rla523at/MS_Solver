#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Matrix.h"
#include "../MS_Solver/INC/Euclidean_Vector.h"
#include "../MS_Solver/INC/Polynomial.h"

TEST(Matrix_Wrapper, change_columns_1)
{
	std::vector<double> value = { 1,2,3,4,5,6 };
	Matrix_Wrapper mw(2, 2, value.data() + 2);
	mw.change_columns(0, 2, 0.0);
	std::vector<double> ref = { 1,2,0,0,0,0 };
	EXPECT_EQ(value, ref);
}
TEST(Matrix_Wrapper, change_columns_2)
{
	std::vector<double> value = { 1,2,3,4,5,6 };
	Matrix_Wrapper mw(2, 3, value.data());
	mw.change_columns(1, 3, 0.0);
	std::vector<double> ref = { 1,2,0,0,0,0 };
	EXPECT_EQ(value, ref);
}

TEST(Matrix, construct_diagonal_matrix1) 
{
	Matrix m(2, {3, 3});
	Matrix ref(2, 2, { 3,0,0,3 });

	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_row1) 
{
	Matrix m(2, 2, { 1,2,3,4 });

	std::vector<double> v = { 5,6 };
	m.change_row(0, v);

	Matrix ref(2, 2, {5, 6, 3, 4});
	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_row2) 
{
	Matrix m(2, 2, { 1,2,3,4 });

	std::vector<double> v = { 5,6 };
	m.change_row(0, v);

	Matrix ref(2, 2, { 5, 6, 3, 4 });
	EXPECT_EQ(m, ref);
}
TEST(Matrix, change_rows1) 
{
	const Matrix m = { 2, 2,  {1,2,3,4} };

	Matrix dm(4, 2);
	dm.change_rows(0, m);

	const Matrix ref(4, 2, { 1,2,3,4,0,0,0,0 });
	EXPECT_EQ(dm, ref);
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
TEST(Matrix, change_columns_3)
{
	Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
	m.change_columns(0, 1, 4);

	const Matrix ref = { 2,3,{ 4,2,3,4,5,6 } };
	EXPECT_EQ(m, ref);
}
TEST(Matrix, column_1)
{
	const Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
	const auto result = m.column(0);

	const std::vector<double> ref = { 1,4 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, column_2)
{
	const Matrix m = { 2,3,{ 1,2,3,4,5,6 } };
	const auto result = m.column(2);

	const std::vector<double> ref = { 3,6 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, get_inverse_1)
{
	Matrix m(2, 2, { 2,1,1,1 });

	const auto result = m.get_inverse();

	Matrix ref(2, 2, { 1,-1,-1,2 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, inverse_1) 
{
	Matrix m(2, 2, { 1,2,3,4 });
	m.inverse();

	Matrix ref(2, 2, { -2,1,1.5,-0.5 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(m.at(i, j), ref.at(i, j));
}
TEST(Matrix, inverse_2)
{
	Matrix m(2, 2, { 1,2,3,4 });
	m.transpose();
	EXPECT_ANY_THROW(m.inverse());

	//Matrix ref(2, 2, { -2,1.5,1,-0.5 });
	//const auto [rows, cols] = ref.size();
	//for (size_t i = 0; i < rows; ++i)
	//	for (size_t j = 0; j < cols; ++j)
	//		EXPECT_DOUBLE_EQ(m.at(i, j), ref.at(i, j));
}
TEST(Matrix, inverse_3) 
{
	Matrix m(5, 5, { 1, 2, 3, 4, 5,  2, 13, 4, 11, 6, 1, 6, 5, 4, 3, 4, 1, 12, 7, 8, 9, 8, 7, 6, 5 });
	m.inverse();

	Matrix ref(5, 5, { -0.03333333333333333, 3.426508252053199e-18, -0.1666666666666667, 2.669049661773008e-17,   0.1333333333333333,              0.09375,               -0.0625,                0.25,                -0.125,              0.03125,  -0.1541666666666667,  -0.04166666666666667,  0.1666666666666667,   0.08333333333333333, -0.02916666666666667,  -0.3395833333333333,    0.2708333333333333, -0.4166666666666667,    0.2083333333333333,             -0.06875,   0.5333333333333333,   -0.1666666666666667,  0.1666666666666667,   -0.1666666666666667,  0.03333333333333333 });
	const auto [rows, cols] = ref.size();
	const double epsilon = 1.0E-15;
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_NEAR(m.at(i, j), ref.at(i, j), epsilon); // suffering by extreme machine error!
}
TEST(Matrix, inverse_4) 
{
	Matrix m(5, 5, { 1, 2, 3, 4, 5,  2, 3, 4, 1, 6, 1, 6, 5, 4, 3, 4, 1, 2, 7, 8, 9, 8, 7, 6, 5 });
	m.inverse();

	Matrix ref(5, 5, { -0.033333333333333,  -0.000000000000000,  -0.166666666666667,  -0.000000000000000,   0.133333333333333,  -1.291666666666667,   0.250000000000000,   0.666666666666667,   0.500000000000000,  -0.208333333333333,   1.616666666666666,  -0.250000000000000,  -0.666666666666667,  -0.750000000000000,   0.283333333333333,   0.275000000000000,  -0.250000000000000,   0.000000000000000,   0.000000000000000,   0.025000000000000,  -0.466666666666667,   0.250000000000000,   0.166666666666667,   0.250000000000000,  -0.133333333333333 });
	const auto [rows, cols] = ref.size();
	const double epsilon = 1.0E-15;
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_NEAR(m.at(i, j), ref.at(i, j), epsilon); // suffering by extreme machine error!
}
TEST(Matrix, inverse_5)
{
	Matrix m(2, 2, { 2,1,1,1 });
	m.inverse();

	Matrix ref(2, 2, { 1,-1,-1,2 });
	EXPECT_EQ(m, ref);
}
TEST(Matrix, operator_addition1)
{
	const Matrix m1 = { 1,2, { 1,2 } };
	const Matrix m2 = { 1,2, { 2,2 } };
	const auto result = m1 + m2;

	const Matrix ref = { 1,2,{ 3,4 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition2) 
{
	const Matrix m1 = { 2,1, { 1,2 } };
	const Matrix m2 = { 2,1, { 2,2 } };
	const auto result = m1 + m2;

	const Matrix ref = { 2,1,{ 3,4 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition3) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	const auto result = m1 + m2;

	Matrix ref(2, 2, { 2,5,7,9 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition4) 
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(3, 2, { 1,3,4,5,5,6 });
	EXPECT_ANY_THROW(m1 + m2);
}
TEST(Matrix, operator_addition5) 
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(3, 2, { 1,3,4,5,5,6 });
	m2.transpose();
	EXPECT_ANY_THROW(m1 + m2);

	//const auto result = m1 + m2;

	//Matrix ref(2, 3, { 2,6,8,7,10,12 });
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition6) 
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(3, 2, { 1,3,4,5,5,6 });
	m2.transpose();
	EXPECT_ANY_THROW(m2 + m1);

	//const auto result = m2 + m1;

	//Matrix ref(2, 3, { 2,6,8,7,10,12 });
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition7) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m1.transpose();
	m2.transpose();
	EXPECT_ANY_THROW(m1 + m2);

	//const auto result = m1 + m2;

	//Matrix ref(2, 2, { 2,7,5,9 });
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition8) 
{
	Matrix m1(2, 5, { 1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345 });
	Matrix m2(5, 2, { 1.234234, 2.3462345, 345.324, 2.6345345, 634523.5, 2345345.3,	 23453.345, 234534.6,	 234523.5, 623452.1 });
	m2.transpose();
	EXPECT_ANY_THROW(m1 + m2);
	//const auto result = m1 + m2;

	//Matrix ref(2, 5, { 2.468734,  347.670345, 634529.8262345,    23461.9124567, 234529.73452345,4.8808545, 9.066986845, 2345347.645345, 234535.834563245,  623459.4245345 });
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_addition9) 
{
	Matrix m1(3, 5, { 1.2345, 2.346345, 6.3262345, 8.5674567, 6.23452345, 2.53462, 6.432452345, 2.345345, 1.234563245, 7.3245345,  789.45978, 74.5789123, 74.23541, 4.7894113, 7894.5134 });
	Matrix m2(5, 3, { 1.234234, 2.3462345, 789456.0,   345.324, 2.6345345, 74.48651,  634523.5, 2345345.3, 710.1846, 23453.345,  234534.6,  12.5487,  234523.5,  623452.1, 421.7456 });
	m2.transpose();
	EXPECT_ANY_THROW(m1 + m2);
	//const auto result = m1 + m2;

	//Matrix ref(3, 5, { 2.468734,  347.670345, 634529.8262345,    23461.9124567, 234529.73452345,    4.8808545, 9.066986845, 2345347.645345, 234535.834563245,  623459.4245345, 790245.45978, 149.0654223,      784.42001,       17.3381113,        8316.259 });
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mc_1)
{
	std::vector<double> ar(9, 1);
	const Matrix m(3, 3, std::move(ar));
	const auto result = m * 10;

	std::vector<double> ref_val(9, 10);

	const Matrix ref(3, 3, std::move(ref_val));
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mc_2)
{
	std::vector<double> ar(100, 1);
	const Matrix m(10, 10, std::move(ar));
	const auto result = m * 10;

	std::vector<double> ref_val(100, 10);

	const Matrix ref(10, 10, std::move(ref_val));
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_1) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 9,13,19,29 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_2) 
{
	Matrix m1(2, 8, { 1,2,3,4,5,6,7,8,8,7,6,5,4,3,2,1 });
	Matrix m2(8, 2, { 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8 });
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 204,204,120,120 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_3) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 13,18,18,26 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_4)
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 7,14,15,32 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_5)
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m1.transpose();
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 10,19,14,28 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_6)
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(2, 3, { 1,2,3,3,2,1 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 14,10,32,28 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_7)
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(2, 3, { 1,2,3,3,2,1 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(3, 3, { 13,10,7,17,14,11,21,18,15 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_8)
{
	Matrix m1(2, 3, { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 });
	Matrix m2(2, 3, { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(3, 3, { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_9)
{
	Matrix m1(2, 3, { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 });
	Matrix m2(2, 3, { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
}
TEST(Matrix, operator_mm_10) 
{
	const Matrix m1 = { 1,2, {1,2} };
	const Matrix m2 = { 2,2,{ 2,2,2,2 } };
	const auto result = m1 * m2;

	const Matrix ref = { 1,2,{ 6,6 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_11) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 9,13,19,29 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_12) 
{
	Matrix m1(2, 8, { 1,2,3,4,5,6,7,8,8,7,6,5,4,3,2,1 });
	Matrix m2(8, 2, { 1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8 });
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 204,204,120,120 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_13) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 13,18,18,26 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_14) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 7,14,15,32 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_15) 
{
	Matrix m1(2, 2, { 1,2,3,4 });
	Matrix m2(2, 2, { 1,3,4,5 });
	m1.transpose();
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 10,19,14,28 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_16) 
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(2, 3, { 1,2,3,3,2,1 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 14,10,32,28 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_17) 
{
	Matrix m1(2, 3, { 1,2,3,4,5,6 });
	Matrix m2(2, 3, { 1,2,3,3,2,1 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(3, 3, { 13,10,7,17,14,11,21,18,15 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_18) 
{
	Matrix m1(2, 3, { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 });
	Matrix m2(2, 3, { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 });
	m1.transpose();
	const auto result = m1 * m2;

	Matrix ref(3, 3, { 18.470224531309999,  16.712738084116999,  10.418514118810000,  22.167911283120002,  21.030398669553001,  14.150019875622000,  27.174477241769999,  26.434492089979003,  18.454029295989997 });
	const auto [rows, cols] = ref.size();
	for (size_t i = 0; i < rows; ++i)
		for (size_t j = 0; j < cols; ++j)
			EXPECT_DOUBLE_EQ(result.at(i, j), ref.at(i, j));
	//EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_19) 
{
	Matrix m1(2, 3, { 1.5479,2.4567123,3.414878,4.41487,5.121,6.15789 });
	Matrix m2(2, 3, { 1.1244,2.48711,3.12314,3.789413,2.9135491,1.264863 });
	m2.transpose();
	const auto result = m1 * m2;

	Matrix ref(2, 2, { 18.515714565372999,  17.342737125037932, 36.932522712599997 , 39.438937931479998 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_20)
{
	const Matrix m1 = { 3,{1,0,0} };
	const Matrix m2 = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const auto result = m1 * m2;

	const Matrix ref = { 3,3,{ 1,2,3,0,0,0,0,0,0 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mm_21)
{
	const Matrix m1 = { 3,{1,0,0} };
	const Matrix m2 = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const auto result = m2 * m1;
	const Matrix ref = { 3,3,{ 1,0,0,4,0,0,7,0,0 } };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mv_1)
{
	const Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	const Euclidean_Vector v = { -1,-1,1 };
	const auto result = m * v;
	const Euclidean_Vector ref = { 0,-3,-6 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mv_2) 
{
	const Matrix m = { 2,2, { 1,2,3,4 } };
	const Euclidean_Vector v = { 1,2 };
	const auto result = m * v;

	const Euclidean_Vector ref = { 5,11 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_mv_3) 
{
	const Matrix m = { 1,2,{ 1,1 } };
	const Euclidean_Vector v = { 1,2 };

	const auto result = (m + m) * v;

	const Euclidean_Vector ref = { 6 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, operator_substitution_1)
{
	Matrix m = { 1,2,{ 1,1 } };
	Matrix result;

	result = std::move(m);
	result.at(0, 1) = 2;

	const Matrix ref = { 1,2,{1,2} };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, row1) 
{
	Matrix dm(2, 3, { 1,2,3,4,5,6 });
	const auto result = dm.row(0);

	const std::vector<double> ref = { 1,2,3 };
	EXPECT_EQ(result, ref);
}
TEST(Matrix, transpose_1)
{
	Matrix result(3, 2, { 1,2,3,4,5,6 });
	result.transpose();

	Matrix ref(2, 3, { 1,3,5,2,4,6 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, transpose_2)
{
	Matrix result(3, 2, { 1,2,3,4,5,6 });
	result.transpose();
	result.transpose();

	Matrix ref(3, 2, { 1,2,3,4,5,6 });
	EXPECT_EQ(result, ref);
}
TEST(Matrix, scalar_multiplication_at_columns_1)
{
	Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	m.scalar_multiplcation_at_columns(0, 1, 3);

	const Matrix ref = { 3,3,{3,2,3,12,5,6,21,8,9} };
	EXPECT_EQ(m, ref);
}
TEST(Matrix, scalar_multiplication_at_columns_2)
{
	Matrix m = { 3,3,{ 1,2,3,4,5,6,7,8,9 } };
	m.scalar_multiplcation_at_columns(1, 3, 3);

	const Matrix ref = { 3,3,{1,6,9,4,15,18,7,24,27} };
	EXPECT_EQ(m, ref);
}








TEST(Matrix_Function, change_column_1)
{
	Polynomial x("x0");
	Polynomial y("x1");

	Vector_Function<Polynomial> vf = { x + y , 2 * x + y };

	Matrix_Function<Polynomial> result(2, 2);
	result.change_column(0, vf[0].gradient());
	result.change_column(1, vf[1].gradient());

	Matrix_Function<Polynomial> ref(2, 2, { 1,2,1,1 });
	EXPECT_EQ(result, ref);
}
//TEST(ms, Jacobian_1) 
//{
//	constexpr ushort domain_dimension = 2;
//
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Vector_Function<Polynomial, domain_dimension> vf = { x * y , 2 * x + y };
//
//	const auto result = ms::Jacobian(vf);
//
//	Matrix_Function<Polynomial, domain_dimension, domain_dimension> ref = { y,x,2,1 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Vector_Function, mv_1) 
//{
//	constexpr ushort domain_dimension = 2;
//
//	Polynomial x("x0");
//	Polynomial y("x1");
//
//	Static_Matrix<2, 2> m = { 1,2,3,4 };
//	Vector_Function<Polynomial<2>, 2> vf = { x , y };
//	const auto result = m * vf;
//
//	Vector_Function<Polynomial<2>, 2> ref = { x + 2 * y ,3 * x + 4 * y };
//	EXPECT_EQ(result, ref);
//}

