#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Euclidean_Vector.h"

TEST(Euclidean_Vector, constructor_1)
{
	std::vector<double> val = { 1,2,3 };
	Euclidean_Vector result = { val.data(), val.data() + val.size() };

	Euclidean_Vector ref = { 1,2,3 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_1) 
{
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 5,7,9 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_2) 
{
	const Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = { 4,5,6,7 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 5,7,9,11 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_3)
{
	const Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = v1;
	const Euclidean_Vector v3 = { 4,5,6,7 };
	const auto result = v1 + v2 + v3;

	const Euclidean_Vector ref = { 6,9,12,15 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_4)
{
	Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = v1;
	const Euclidean_Vector v3 = { 4,5,6,7 };

	v1 *= 0;
	const auto result = v2 + v3;

	const Euclidean_Vector ref = { 5,7,9,11 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_5)
{
	Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = std::move(v1);
	const Euclidean_Vector v3 = { 4,5,6,7 };

	v1 = { 2,3,4,5 };
	const auto result = v2 + v3;

	const Euclidean_Vector ref = { 5,7,9,11 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_addition_assign_1)
{
	Euclidean_Vector v1 = { 1,2,3,4 };
	const Euclidean_Vector v2 = { 4,5,6,7 };
	v1 += v2;

	const Euclidean_Vector ref = { 5,7,9,11 };
	EXPECT_EQ(v1, ref);
}
TEST(Euclidean_Vector, operator_substraction_1) {
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 - v2;

	const Euclidean_Vector ref = { -3,-3,-3 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_scalar_multiplication_1) 
{
	const Euclidean_Vector v1 = { 1,2,3 };
	const auto result = v1 * 2;

	const Euclidean_Vector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, operator_scalar_multiplication_2) 
{
	const Euclidean_Vector v1 = { 1,2,3 };
	const auto result = 2 * v1;

	const Euclidean_Vector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, inner_product_1) 
{
	const Euclidean_Vector v1 = { 1,2,3 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1.inner_product(v2);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, inner_product_2)
{
	Euclidean_Vector v1 = { 3,4 };
	const auto result = v1.inner_product(v1);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, L2_norm_1) 
{
	const Euclidean_Vector v1 = { 3,4 };
	const auto result = v1.L2_norm();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector, normalize_1)
{
	Euclidean_Vector v1 = { 3,4 };
	v1.normalize();

	Euclidean_Vector ref = { 3.0 / 5.0, 4.0 / 5.0 };

	for (int i = 0; i < v1.size(); ++i)
		EXPECT_DOUBLE_EQ(v1[i], ref[i]);
}

TEST(Euclidean_Vector_Wrapper, operator_addition_1)
{
	const std::vector<double> vec = { 1,2,3 };
	const Euclidean_Vector_Wrapper v1 = vec;
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 5,7,9 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, operator_addition_2)
{
	std::vector<double> vec = { 1,2,3 };
	const Euclidean_Vector_Wrapper v1 = vec;
	vec = { 4,3,2 };
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1 + v2;

	const Euclidean_Vector ref = { 8,8,8 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, operator_scalar_multiplication_1)
{
	const std::vector<double> vec = { 1,2,3 };
	const Euclidean_Vector_Wrapper v1 = vec;
	const auto result = 2 * v1;

	const Euclidean_Vector ref = { 2,4,6 };
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, inner_product_1)
{
	const std::vector<double> vec = { 1,2,3 };
	const Euclidean_Vector_Wrapper v1 = vec;
	const Euclidean_Vector v2 = { 4,5,6 };
	const auto result = v1.inner_product(v2);

	const double ref = 32;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, inner_product_2)
{
	std::vector<double> vec = { 3,4 };
	const Euclidean_Vector_Wrapper v1 = vec;
	const auto result = v1.inner_product(v1);

	const auto ref = 25;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, L1_norm_1)
{
	std::vector<double> vec = { -3,4 };
	const Euclidean_Vector_Wrapper v1 = vec;

	const auto result = v1.L1_norm();

	const auto ref = 7;
	EXPECT_EQ(result, ref);
}
TEST(Euclidean_Vector_Wrapper, L2_norm_1)
{
	std::vector<double> vec = { 3,4 };
	const Euclidean_Vector_Wrapper v1 = vec;

	const auto result = v1.L2_norm();

	const auto ref = 5;
	EXPECT_EQ(result, ref);
}

//TEST(Dynamic_Euclidean_Vector, be_absolute_1) 
//{
//	Dynamic_Euclidean_Vector v1 = { -1, -2 ,3 ,-4 ,5 };
//	v1.be_absolute();
//
//	const Dynamic_Euclidean_Vector ref = { 1,2,3,4,5 };
//	EXPECT_EQ(v1, ref);
//}

//
//TEST(ms, min_value_gathering_vector_1) 
//{
//	Euclidean_Vector v1 = { 1,2,3 };
//	Euclidean_Vector v2 = { 2,3,1 };
//	Euclidean_Vector v3 = { 3,1,2 };
//
//	std::vector<Euclidean_Vector<3>> vec = { v1,v2,v3 };
//
//	const auto result = ms::gather_min_value(vec);
//	const Euclidean_Vector ref = { 1,1,1 };
//	EXPECT_EQ(result, ref);
//}
//
//TEST(ms, max_value_gathering_vector_1) 
//{
//	Euclidean_Vector v1 = { 1,2,3 };
//	Euclidean_Vector v2 = { 2,3,1 };
//	Euclidean_Vector v3 = { 3,2,1 };
//
//	std::vector<Euclidean_Vector<3>> vec = { v1,v2,v3 };
//
//	const auto result = ms::gather_max_value(vec);
//	const Euclidean_Vector ref = { 3,3,3 };
//	EXPECT_EQ(result, ref);
//}