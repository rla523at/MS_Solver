#include "gtest/gtest.h"
#include "../MS_Solver/INC/Vector_Function.h"
#include "../MS_Solver/INC/Polynomial.h"

TEST(Vector_Function, constructor) {	
	constexpr ushort domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	constexpr ushort range_dimension = 3;
	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf = { 1,x,y };

	EXPECT_EQ(vf.at(0), 1);
	EXPECT_EQ(vf.at(1), x);
	EXPECT_EQ(vf.at(2), y);
}

TEST(Vector_Function, mv_1) {
	constexpr ushort domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Matrix<2, 2> m = { 1,2,3,4 };
	Vector_Function<Polynomial<2>, 2> vf = { x , y };
	const auto result = m * vf;

	Vector_Function<Polynomial<2>, 2> ref = { x + 2 * y ,3 * x + 4 * y };
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Vector_Function, differentiate_1)
{
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>,domain_dimension> vf = { x + y , 2 * x + y };
	const auto result = vf.differentiate(0);

	Vector_Function<Polynomial<domain_dimension>, domain_dimension> ref = { 1,2 };
	EXPECT_EQ(ref, result);
}
GTEST_TEST(Vector_Function, differentiate_2){
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>, domain_dimension> vf = { 0.25 * x * y + 1.25 * x + 0.25 * y + 2.25, -0.25 * x * y - 0.75 * x + 0.25 * y + 1.75 };
	const auto result = vf.differentiate(0);

	Vector_Function<Polynomial<domain_dimension>, domain_dimension> ref = { 0.25 * y + 1.25,-0.25 * y - 0.75 };
	EXPECT_EQ(ref, result);
}

GTEST_TEST(Vector_Function, cross_product_1){
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	constexpr ushort range_dimension = 3;

	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf1 = { 1 + x ,		2 * x + y,		3 };
	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf2 = { 2 * x - y ,	1 + (-1 * x),	2 };
	const auto result = vf1.cross_product(vf2);

	Vector_Function<Polynomial<domain_dimension>, range_dimension> ref = { 7 * x + 2 * y - 3,4 * x - 3 * y - 2, 1 - 5 * (x ^ 2) + (y ^ 2) };
	EXPECT_EQ(ref, result);
}
GTEST_TEST(Vector_Function, cross_product_2) {
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	constexpr ushort range_dimension = 2;

	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf1 = { 1 + x ,		2 * x + y };
	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf2 = { 2 * x - y ,	1 + (-1 * x) };
	const auto result = vf1.cross_product(vf2);

	Vector_Function<Polynomial<domain_dimension>, 3> ref = { 0, 0, 1 - 5 * (x ^ 2) + (y ^ 2) };
	EXPECT_EQ(ref, result);
}

TEST(Vector_Function, L2_norm_1) {
	constexpr size_t domain_dimension = 2;
	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	constexpr ushort range_dimension = 3;

	Vector_Function<Polynomial<domain_dimension>, range_dimension> vf1 = { 1, x, y };
	const auto result = vf1.L2_norm();

	Irrational_Function<domain_dimension> ref(1 + (x ^ 2) + (y ^ 2), 0.5);
	EXPECT_EQ(ref, result);
}


GTEST_TEST(Dynamic_Vector_Function, at_1)
{
	Polynomial<2> x("x0");
	Polynomial<2> y("x1");

	Dynamic_Vector_Function<Polynomial<2>> vf = { 0 , 2 * x + y };

	const auto ref = 2 * x + y;
	EXPECT_EQ(ref, vf.at(1));
}

GTEST_TEST(Dynamic_Vector_Function, operator_call_1)
{
	constexpr size_t domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Dynamic_Vector_Function<Polynomial<domain_dimension>> vf = { x + y , 2 * x + y };
	Euclidean_Vector v = { 1,1 };
	const auto result = vf(v);

	Dynamic_Euclidean_Vector ref = { 2,3 };
	EXPECT_EQ(ref, result);
}

TEST(Matrix_Function, change_column) {
	constexpr ushort domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>, domain_dimension> vf = { x + y , 2 * x + y };
		
	Matrix_Function<Polynomial<domain_dimension>, domain_dimension, domain_dimension> result;
	result.change_column(0, vf[0].gradient());
	result.change_column(1, vf[1].gradient());

	Matrix_Function<Polynomial<domain_dimension>, domain_dimension, domain_dimension> ref = { 1,2,1,1 };
	EXPECT_EQ(result, ref);
}

TEST(ms, Jacobian_1) {
	constexpr ushort domain_dimension = 2;

	Polynomial<domain_dimension> x("x0");
	Polynomial<domain_dimension> y("x1");

	Vector_Function<Polynomial<domain_dimension>, domain_dimension> vf = { x * y , 2 * x + y };

	const auto result = ms::Jacobian(vf);

	Matrix_Function<Polynomial<domain_dimension>, domain_dimension, domain_dimension> ref = { y,x,2,1 };
	EXPECT_EQ(result, ref);

}