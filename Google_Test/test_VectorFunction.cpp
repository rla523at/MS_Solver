
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Vector_Function.h"
#include "../MS_Solver/INC/Polynomial.h"

TEST(Vector_Function, OPERATOR_CALL1) 
{
	Polynomial x("x0");
	Polynomial y("x1");
	Polynomial z("x2");

	auto p1 = (x ^ 2) + 3 * (x ^ 2) * y + (y ^ 3) + (z ^ 2) - 6;
	auto p2 = x + y + z - 3;
	auto p3 = (y ^ 2) * z + x * z - 2;
	
	Vector_Function<Polynomial> f = { p1,p2,p3 };
	std::vector<double> node = { 1,1,1 };
	const auto result = f(node);

	std::vector<double> ref = { 0,0,0 };
	EXPECT_EQ(result, ref);
}
//TEST(Vector_Function, DIFFERENTIATE1) {
//	Polynomial x("x0");
//	Polynomial y("x1");
//	Polynomial z("x2");
//
//	auto p1 = x + 1;
//	auto p2 = x + y + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 1,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(Vector_Function, DIFFERENTIATE2) {
//	auto p1 = x * y + 1;
//	auto p2 = x + y * z + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { y,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(Vector_Function, DIFFERENTIATE3) {
//	auto p1 = x * y + 1;
//	auto p2 = x + y * z + z - 3;
//	auto p3 = y + z;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { x,z,1 };
//	EXPECT_EQ(f, ref);
//}
//TEST(Vector_Function, DIFFERENTIATE4) {
//	auto p1 = 1.5 * x + 0.5 * y + 3;
//	auto p2 = y + 3;
//	auto p3 = 0;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 1.5,0,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(Vector_Function, DIFFERENTIATE5) {
//	auto p1 = 1.5 * x + 0.5 * y + 3;
//	auto p2 = y + 3;
//	auto p3 = 0;
//	Vector_Function<Polynomial> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	Vector_Function<Polynomial> ref = { 0.5,1,0 };
//	EXPECT_EQ(f, ref);
//}
TEST(Vector_Function, CROSS_PRODUCT1) 
{
	Vector_Function<Polynomial> vf1 = { 1.5,0,0 };
	Vector_Function<Polynomial> vf2 = { 0.5,1,0 };
	const auto result = vf1.cross_product(vf2);

	Vector_Function<Polynomial> ref = { 0,0,1.5 };
	EXPECT_EQ(result, ref);
}
//TEST(Vector_Function, L2_NORM1) 
//{
//	Vector_Function<Polynomial> vf = { 0,0,1.5 };
//	const auto result = vf.L2_norm();
//
//	Polynomial ref = 1.5;
//	EXPECT_EQ(result, ref);
//}
