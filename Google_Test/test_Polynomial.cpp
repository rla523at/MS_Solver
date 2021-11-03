#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Polynomial.h"

TEST(Simple_Poly_Term, operator_addition_1) {
	const auto spt1 = Simple_Poly_Term(0.0);
	const auto spt2 = Simple_Poly_Term("x0");
	const auto result = spt1 + spt2;

	const auto ref = Simple_Poly_Term("x0");
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_addition_2) {
	const auto spt1 = Simple_Poly_Term(0.0);
	const auto spt2 = Simple_Poly_Term("x3");
	const auto result = spt1 + spt2;

	std::vector<double> coefficients = { 0,0,0,1 };
	const auto ref = Simple_Poly_Term(coefficients);
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_addition_3) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, 1);

	std::vector<double> coefficients2 = { 3,2,1 };
	const auto spt2 = Simple_Poly_Term(coefficients2, -3);

	const auto result = spt1 + spt2;

	std::vector<double> coefficients_ref = { 4,4,4 };
	const auto ref = Simple_Poly_Term(coefficients_ref, -2);

	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_addition_4) {
	const auto spt1 = Simple_Poly_Term({ 1 }, -4);
	const auto spt2 = Simple_Poly_Term({ -1 }, 4);
	const auto result = spt1 + spt2;

	const Simple_Poly_Term ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_call_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, 1);

	std::vector<double> domain_vector = { 1,2,3 };
	const auto result = spt(domain_vector);

	const auto ref = 15;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_call_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, 1);

	std::array<double,3> domain_vector = { 1,2,3 };
	const auto result = spt(domain_vector);

	const auto ref = 15;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_call_3) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);

	std::vector<double> domain_vector = { 1,2,3,4,5,6 };
	const auto result = spt(domain_vector);

	const auto ref = 10;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, differentiate_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);

	constexpr ushort variable_index = 0;
	const auto result = spt.differentiate(variable_index);

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, differentiate_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);

	constexpr ushort variable_index = 2;
	const auto result = spt.differentiate(variable_index);

	const auto ref = 3;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, differentiate_3) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);

	constexpr ushort variable_index = 5;
	const auto result = spt.differentiate(variable_index);

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Simple_Poly_Term, operator_less_then_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, -4);

	std::vector<double> coefficients2 = { 1,2,3 };
	const auto spt2 = Simple_Poly_Term(coefficients2, -3);

	EXPECT_TRUE(spt1 < spt2);
}
TEST(Simple_Poly_Term, operator_less_then_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, -4);

	std::vector<double> coefficients2 = { 1,2,3 };
	const auto spt2 = Simple_Poly_Term(coefficients2, -5);

	EXPECT_FALSE(spt1 < spt2);
}
TEST(Simple_Poly_Term, operator_less_then_3) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, -4);

	std::vector<double> coefficients2 = { 1,2,3 };
	const auto spt2 = Simple_Poly_Term(coefficients2, -4);

	EXPECT_FALSE(spt1 < spt2);
}
TEST(Simple_Poly_Term, operator_less_then_4) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, -4);
	const auto spt2 = Simple_Poly_Term("x3");

	EXPECT_TRUE(spt1 < spt2);
}
TEST(Simple_Poly_Term, operator_less_then_5) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt1 = Simple_Poly_Term(coefficients1, -4);
	
	std::vector<double> coefficients2 = { 1,3,2 };
	const auto spt2 = Simple_Poly_Term(coefficients2, -4);

	EXPECT_TRUE(spt1 < spt2);
}
TEST(Simple_Poly_Term, domain_dimension_1) {	
	const auto spt1 = Simple_Poly_Term({ 1 }, -4);	
	const auto spt2 = Simple_Poly_Term({ -1 }, 4);
	const auto spt3 = spt1 + spt2;
	const auto result = spt3.domain_dimension();

	const auto ref = 0;
	EXPECT_EQ(result,ref);
}
TEST(Simple_Poly_Term, domain_dimension_2) {
	const auto y = Simple_Poly_Term("x1");
	const auto result = y.domain_dimension();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}

TEST(PoweredPolyTerm, operator_scalar_multimplication_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);
	const auto ppt = PoweredPolyTerm(spt, 2);
		
	constexpr auto scalar = 3;
	const auto result = ppt * scalar;

	const auto ref = PolyTerm(scalar,ppt);
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, operator_scalar_multimplication_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);
	const auto ppt = PoweredPolyTerm(spt, 2);

	constexpr auto scalar = -1.908;
	const auto result = ppt * scalar;

	const auto ref = PolyTerm(scalar, ppt);
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, operator_call_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, -4);
	const auto ppt = PoweredPolyTerm(spt, 2);

	std::vector<double> domain_vector = { 1,2,3,4,5,6 };
	const auto result = ppt(domain_vector);

	const auto ref = 100;
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, operator_call_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, 1);
	const auto ppt = PoweredPolyTerm(spt, 3);

	std::vector<double> domain_vector = { 1,2,3 };
	const auto result = ppt(domain_vector);

	constexpr auto ref = 3375;
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, differentiate_1) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, 1);
	const auto ppt = PoweredPolyTerm(spt, 3);

	constexpr ushort variable_index = 0;
	const auto result = ppt.differentiate(variable_index);
	
	const auto ppt2 = PoweredPolyTerm(spt, 2);
	const auto ref = PolyTerm(3, ppt2);
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, differentiate_2) {
	std::vector<double> coefficients1 = { 1,2,3 };
	const auto spt = Simple_Poly_Term(coefficients1, 1);
	const auto ppt = PoweredPolyTerm(spt, 3);

	constexpr ushort variable_index = 2;
	const auto result = ppt.differentiate(variable_index);

	const auto ppt2 = PoweredPolyTerm(spt, 2);
	const auto ref = PolyTerm(9, ppt2);
	EXPECT_EQ(result, ref);
}
TEST(PoweredPolyTerm, differentiate_3) {
	const auto spt1 = Simple_Poly_Term("x0");
	const auto spt2 = spt1 + 1;

	const auto ppt1 = PoweredPolyTerm(spt2, 2);

	constexpr auto variable_index = 0;
	const auto result = ppt1.differentiate(variable_index);

	const auto ref = PolyTerm(2, spt2);
	EXPECT_EQ(result, ref);
}

TEST(PolyTerm, operator_multiplication_1) {
	std::vector<double> coefficients = { 1 };
	const auto spt = Simple_Poly_Term(coefficients, 1);
	const auto ppt = PoweredPolyTerm(spt, 2);

	const auto pt1 = PolyTerm(3, ppt);
	const auto pt2 = PolyTerm(4, ppt);
	const auto result = pt1 * pt2;

	const auto ref_ppt = PoweredPolyTerm(spt, 4);
	const auto ref = PolyTerm(12, ref_ppt);
	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, operator_call_1) {
	std::vector<double> coefficients = { 1 };
	const auto spt1 = Simple_Poly_Term(coefficients, 1);
	const auto spt2 = Simple_Poly_Term(coefficients, 2);
	const auto spt3 = Simple_Poly_Term(coefficients, 3);
	const auto spt4 = Simple_Poly_Term(coefficients, 4);
	const auto spt5 = Simple_Poly_Term(coefficients, 5);

	const auto pt1 = PolyTerm(spt1);
	const auto pt2 = PolyTerm(spt2);
	const auto pt3 = PolyTerm(spt3);
	const auto pt4 = PolyTerm(spt4);
	const auto pt5 = PolyTerm(spt5);

	const auto pt = pt1 * pt2 * pt3 * pt4 * pt5;

	std::vector<double> domain_vector = { 0 };
	const auto result = pt(domain_vector);

	const auto ref = 120;
	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, operator_call_2) {
	std::vector<double> coefficients = { 1 };
	const auto spt1 = Simple_Poly_Term(coefficients, 1);
	const auto spt2 = Simple_Poly_Term(coefficients, 2);
	const auto spt3 = Simple_Poly_Term(coefficients, 3);
	const auto spt4 = Simple_Poly_Term(coefficients, 4);
	const auto spt5 = Simple_Poly_Term(coefficients, 5);

	const auto pt1 = PolyTerm(spt1);
	const auto pt = pt1 * spt2 * spt3 * spt4 * spt5;

	std::vector<double> domain_vector = { 0 };
	const auto result = pt(domain_vector);

	const auto ref = 120;
	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, domain_dimension_1) {
	const auto x = Simple_Poly_Term("x0");
	const auto y = Simple_Poly_Term("x1");

	const auto pt1 = x * y;

	const auto result = pt1.domain_dimension();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, domain_dimension_2) {
	const auto x = Simple_Poly_Term("x0");
	const auto y = Simple_Poly_Term("x1");

	const auto pt1 = (x + 1) * x;

	const auto result = pt1.domain_dimension();

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, differentiate_1) {
	const auto spt1 = Simple_Poly_Term("x0");	
	const auto spt2 = spt1 + 1;

	const auto pt1 = PolyTerm(spt1);
	const auto pt2 = PolyTerm(spt2);
	const auto pt3 = pt1 * pt2;

	constexpr size_t variable_index = 0;
	const auto result = pt3.differentiate(variable_index);

	const auto x = Polynomial("x0");
	const auto ref = 2 * x + 1;

	EXPECT_EQ(result, ref);
}
TEST(PolyTerm, differentiate_2) {
	const auto spt1 = Simple_Poly_Term("x0");
	const auto spt2 = spt1 + 1;


	const auto pt1 = PolyTerm(spt2);
	const auto pt2 = pt1 * pt1;

	constexpr size_t variable_index = 0;
	const auto result = pt2.differentiate(variable_index);

	const auto x = Polynomial("x0");
	const auto ref = 2 * x + 2;

	EXPECT_EQ(result, ref);
}

TEST(Polynomial, operator_addition_1) {
	Polynomial x("x0");
	Polynomial y("x1");

	const auto p1 = 3 * x + y;
	const auto p2 = 5 * (x ^ 2);
	const auto p3 = y ^ 2;
	const auto p4 = -7;
	const auto result = p1 + p2 + p3 + p4;

	const auto ref = 5 * (x ^ 2) + (y ^ 2) + 3 * x + y - 7;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_2) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + 1;
	const auto p2 = X + Y + 1;
	const auto result = p1 + p2;

	const auto ref = 2 * X + Y + 2;
	EXPECT_EQ(result, ref);	
}
TEST(Polynomial, operator_addition_3) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X ^ 2) + X;
	const auto p2 = X + Y + 1;
	const auto result = p1 + p2;

	auto ref = (X ^ 2) + 2 * X + Y + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_4) {
	const Polynomial p1 = 0.0;
	const Polynomial p2 = 0.0;
	const auto result = p1 + p2;

	const Polynomial ref = 0.0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_5) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = -12 * (X ^ 2) * Y + 45 * X;
	const auto p2 = 12 * (X ^ 2) * Y - 45 * X;
	const auto result = p1 + p2;

	const Polynomial ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_6) {
	Polynomial X("x0");
	
	const auto p1 = (X ^ 2) + X;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = (X ^ 2) + 2 * X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_7) {
	Polynomial X("x0");

	const Polynomial p1 = 0.0;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_8) {
	Polynomial X("x0");

	const Polynomial p1 = 0.0;
	const auto p2 = X + 1;
	const auto result = p2 + p1;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_9) {
	Polynomial X("x0");

	const Polynomial p1 = 3;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = X + 4;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_10) {
	Polynomial X("x0");
	Polynomial Y("x1");
	Polynomial Z("x2");

	const auto p1 = X + Y + Z + 5;
	const auto p2 = X + 1;
	const auto p3 = X + Y + 3;
	const auto p4 = X + Y + Z - 9;
	const auto result = p1 + p2 + p3 + p4;

	const auto ref = 4 * X + 3 * Y + 2 * Z;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_addition_11) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X;
	const auto p2 = Y;
	const auto result = p1 + p2;

	const auto ref = X + Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_1) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + Y;
	const auto p2 = X + 3 * Y;
	const auto result = p1 - p2;

	const auto ref = -2 * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_2) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X * Y;
	const auto p2 = X + 3 * Y;
	const auto result = p1 - p2;

	auto ref = X * Y - X - 3 * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_3) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X * Y;
	const auto p2 = X;
	const auto result = p1 - p2;

	const auto ref = X * Y - X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_4) {
	Polynomial X("x0");

	const Polynomial p1 = 0;
	const auto p2 = X;
	const auto result = p1 - p2;

	const auto ref = -1 * X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_5) {
	Polynomial X("x0");
	Polynomial Y("x1");
	Polynomial Z("x2");

	const Polynomial p1 = 0;
	const auto p2 = 3 * (X * Y) - 3 * ((X ^ 2) * Z);
	const auto result = p1 - p2;

	const auto ref = -3 * (X * Y) + 3 * ((X ^ 2) * Z);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_6) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X - 1) * (Y - 1);
	const Polynomial p2 = 1;
	const auto result = p1 - p2;

	const auto ref = X * Y - X - Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_7) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = X + 1;
	const auto result = p1 - p2;

	const auto ref = X * Y - 2 * X - Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_8) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = 2 * ((X - 1) ^ 2) + X + Y - 1;
	const auto result = p1 - p2;

	const auto ref = -2 * (X ^ 2) + X * Y + 2 * X - 2 * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_9) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = (X - 1) * (Y - 1);
	const auto result = p1 - p2;

	const auto ref = 0.0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_substraction_10) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = X * Y - Y;
	const auto result = p1 - p2;

	const auto ref = -1 * X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_1) {
	Polynomial X("x0");

	const auto p1 = X + 1;
	const auto p2 = X - 2;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) - X - 2;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_2) {
	Polynomial X("x0");

	const auto p1 = (X ^ 2) - 1;
	const auto p2 = (X ^ 2) + 2 * X + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) + 2 * (X ^ 3) - 2 * X - 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_3) {
	Polynomial X("x0");

	const auto p1 = (X ^ 2) + 2 * X + 1;
	const auto p2 = (X ^ 2) - 2 * X + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) - 2 * (X ^ 2) + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_4) {
	const Polynomial p1 = 3;
	const Polynomial p2 = 2;
	const auto result = p1 * p2;

	const Polynomial ref = 6;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_5) {
	Polynomial X("x0");

	const Polynomial p1 = 3;
	const auto p2 = X - 1;
	const auto result = p1 * p2;

	const auto ref = 3 * X - 3;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_6) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X;
	const auto p2 = X + Y + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + X * Y + X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_7) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + Y + 1;
	const auto p2 = X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + X * Y + X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_8) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + Y;
	const auto p2 = Y + X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + 2 * X * Y + (Y ^ 2);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_9) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + Y;
	const auto result = -1 * p1;

	const auto ref = -1 * X - 1 * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_10) {
	Polynomial X("x0");

	const auto p1 = X + 1;
	const auto p2 = X + 1;
	const auto result = p1 * p2;

	auto ref = (X ^ 2) + 2 * X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_11) {
	const Polynomial p1 = 0;
	const Polynomial p2 = 0;
	const auto result = p1 * p2;

	const Polynomial ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_12) {
	Polynomial X("x0");

	const auto p1 = 1;
	const auto p2 = X + 1;
	const auto result = p1 * p2;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_13) {
	Polynomial X("x0");

	const auto p1 = X + 1;
	const Polynomial p2 = 0;
	const auto result = p1 * p2;

	const Polynomial ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_14) {
	Polynomial X("x0");

	const auto p1 = (X ^ 2) + X + 1;
	const auto result = 5 * p1;

	const auto ref = 5 * (X ^ 2) + 5 * X + 5;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_15) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = 3 * X;
	const auto p2 = 2 * Y;
	const auto result = p1 * p2;

	const auto ref = 6 * X * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_16) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = 3 * X * Y;
	const auto p2 = (2 * X) ^ 2;
	const auto result = p1 * p2;

	const auto ref = 12 * (X ^ 3) * Y;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_17) {
	Polynomial X("x0");
	Polynomial Y("x1");
	Polynomial Z("x2");

	const auto p1 = X * (Y ^ 2) * (Z ^ 3);
	const auto p2 = X * (Y ^ 2);
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) * (Y ^ 4) * (Z ^ 3);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_18) {
	Polynomial X("x0");
	Polynomial Y("x1");
	Polynomial Z("x2");

	const auto p1 = X * (Y ^ 4);
	const auto p2 = X * (Z ^ 3);
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) * (Y ^ 4) * (Z ^ 3);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_19) {
	Polynomial X("x0");
	Polynomial Y("x1");
	
	const auto p1 = X * (X + Y);
	const auto p2 = (Y + X) * X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) + 2 * (X ^ 3) * Y + ((X * Y) ^ 2);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_20) {
	Polynomial X("x0");
	Polynomial Y("x1");
	
	const auto p1 = 0.25*Y +1.25;
	const auto p2 = -0.25 * X + 0.25;
	const auto result = p1 * p2;

	const auto ref = (-1 * X * Y - 5 * X + Y + 5) * (1.0 / 16.0);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_21) {
	Polynomial X("x0");
	Polynomial Y("x1");
	
	const auto p1 = X + 1;
	const auto p2 = Y + 1;
	const auto p3 = X - 1;
	const auto p4 = Y - 1;
	const auto result = p1 * p2 * p3 * p4;

	const auto ref = (X ^ 2) * (Y ^ 2) - 1 * (X ^ 2) - 1 * (Y ^ 2) + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_22) {	
	Polynomial Y("x1");

	const auto p1 = Y + 1;
	const auto p2 = Y - 1;
	const auto result = p1 * p2;

	const auto ref = (Y ^ 2) - 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_multiplication_23) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + 1;
	const auto p2 = Y + 1;
	const auto p3 = X - 1;
	const auto p4 = p1 * p2 * p3;
	const auto result = p4 * (Y - 1);

	const auto ref = (X ^ 2) * (Y ^ 2) - 1 * (X ^ 2) - 1 * (Y ^ 2) + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_call_1) {
	const Polynomial m = 1;
	const auto values = { 1,2,3 };
	const auto result = m(values);

	constexpr double ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_call_2) {
	const Polynomial m = 1;
	const std::vector<double> values = { 0,0,0 };
	const auto result = m(values);

	constexpr double ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, operator_call_3) {
	const auto X = Polynomial("x0");
	const auto Y = Polynomial("x1");
	const auto Z = Polynomial("x2");

	const auto p = 3 * X * (Y ^ 2) * (Z ^ 3);
	const std::vector<double> values = { 1,2,3 };
	const auto result = p(values);

	constexpr double ref2 = 324;
	EXPECT_EQ(result, ref2);
}
TEST(Polynomial, operator_call_4) {
	const auto X = Polynomial("x0");
	const auto Y = Polynomial("x1");
	const auto Z = Polynomial("x2");

	const auto p = X * (Y ^ 2) * (Z ^ 3);
	const std::vector<double> values = { 1.84,2.789,3.487946 };
	const auto result = p(values);

	constexpr double ref3 = 6.073291260986822e+02;
	EXPECT_DOUBLE_EQ(result, ref3);
}
TEST(Polynomial, operator_call_5) {
	constexpr size_t space_dimension = 3;

	const auto x = Polynomial("x0");
	const auto y = Polynomial("x1");
	const auto z = Polynomial("x2");

	auto p = x + y + z;
	p.be_absolute();

	const std::vector<double> values = { 4,-7,-8 };
	const auto result = p(values);

	constexpr auto ref = 11;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Polynomial, operator_call_6) {
	constexpr size_t space_dimension = 3;

	const auto x = Polynomial("x0");
	const auto y = Polynomial("x1");
	const auto z = Polynomial("x2");

	auto p = x + y + z;
	const std::vector<double> values = { 4,-7,-8 };
	const auto result = p(values);

	constexpr auto ref = -11;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Polynomial, operator_power_1) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p = -0.125 * X + 0.125 * Y + 0.5;
	const auto result = p ^ 2;

	const auto ref = 0.015625 * (X ^ 2) - 0.03125 * X * Y + 0.015625 * (Y ^ 2) - 0.125 * X + 0.125 * Y + 0.25;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_1) {
	Polynomial X("x0");

	const auto p1 = X + 1;
	const auto p2 = X + 2;
	const auto p3 = X + 3;
	const auto p4 = p1 * p2 + p3;
	const auto p5 = p4 * p4;
	const auto result = 2 * (p5 + p1);

	auto ref = 2 * (X ^ 4) + 16 * (X ^ 3) + 52 * (X ^ 2) + 82 * X + 52;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_2) {
	Polynomial X("x0");

	auto p1 = X + 1;
	const auto p2 = p1 ^ 2;
	const auto p3 = p2 + X + 2;
	const auto p4 = p3 + 2 * X + 3;
	const auto result = 5 * p4;

	const auto ref = 5 * (X ^ 2) + 25 * X + 30;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_3) {
	Polynomial X("x0");
	Polynomial Y("x1");

	auto p1 = X + Y;
	auto p2 = Y + X;
	auto p3 = X;
	const auto result = (p1 * p3) * (p2 * p3);

	auto ref = (X ^ 4) + 2 * (X ^ 3) * Y + ((X * Y) ^ 2);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_4) {
	Polynomial X("x0");
	Polynomial Y("x1");

	auto p1 = 0.25 * Y + 1.25;
	auto p2 = -0.25 * Y + -0.75;
	auto p3 = 0.25 * X + 0.25;
	auto p4 = -0.25 * X + 0.25;
	const auto result = (p1 * p4) - (p2 * p3);

	auto ref = -0.125 * X + 0.125 * Y + 0.5;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_5) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X + 1;
	const auto p2 = Y + 1;
	const auto p3 = X - 1;
	const auto p4 = Y - 1;
	const auto p5 = p1 * p2 - p3 * p4;
	const auto result = p5 ^ 2;

	const auto ref = 4 * (X ^ 2) + 8 * X * Y + 4 * (Y ^ 2);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, complex_operation_6) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = 0.25 * Y + 1.25;
	const auto p2 = -0.25 * Y - 0.75;
	const auto p3 = 0.25 * X + 0.25;
	const auto p4 = -0.25 * X + 0.25;
	const auto p5 = p1 * p4 - p2 * p3;
	const auto result = p5 ^ 2;

	const auto ref = ((X ^ 2) + (Y ^ 2) - 2 * X * Y - 8 * X + 8 * Y + 16) * (1.0 / 64.0);
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, domain_dimension_1) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = -12 * (X ^ 2) * Y + 45 * X;
	const auto p2 = 12 * (X ^ 2) * Y - 45 * X;
	const auto p3 = p1 + p2;
	const auto result = p3.domain_dimension();

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, domain_dimension_2) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X * Y;
	const auto result = p1.domain_dimension();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, domain_dimension_3) {
	Polynomial X("x0");

	const auto p1 = X + 1;
	const auto p2 = X - 2;
	const auto p3 = p1 * p2;
	const auto result = p3.domain_dimension();

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, order_1) {
	const Polynomial p = 0;
	const auto result = p.degree();

	constexpr size_t ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, order_2) {
	const Polynomial X("x0");

	const auto p = X;
	const auto result = p.degree();

	constexpr size_t ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, order_3) {
	const Polynomial X("x0");
	const Polynomial Y("x1");

	const auto p = (X * Y)^2;
	const auto result = p.degree();

	constexpr size_t ref = 4;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_1) {
	Polynomial X("x0");

	const auto p = (X + 1) * (X - 1) + 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	const auto ref = 2 * X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_2) {
	Polynomial X("x0");

	const auto p = (X - 1) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	const auto ref = 2 * X - 2;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_3) {
	Polynomial X("x0");

	const auto p = (X - 1) ^ 2;

	constexpr size_t variable_index = 1;
	const auto result = p.differentiate(variable_index);

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_4) {
	Polynomial X("x0");

	const auto p1 = (X + 1) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = 2 * X + 2;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_5) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X + Y + 2) * (X + 1);

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = 2 * X + Y + 3;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_6) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X + Y + 2) * (X + 1);

	constexpr size_t variable_index = 1;
	const auto result = p1.differentiate(variable_index);

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_7) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = (X + 1) * (Y + 1) + X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = Y + 2;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_8) {
	Polynomial X("x0");

	const auto p1 = ((X + 1) ^ 2) * (X + 2);

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = 3 * (X ^ 2) + 8 * X + 5;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_9) {
	Polynomial X("x0");

	const auto p = (2 * X + 3) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	const auto ref = 8 * X + 12;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_10) {
	Polynomial X("x0");

	const auto p = ((2 * X + 3) ^ 2) * X;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	const auto ref = 12 * (X ^ 2) + 24 * X + 9;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_11) {
	Polynomial X("x0");

	const auto  p = X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	Polynomial ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_12) {
	Polynomial X("x0");

	const auto p = (X ^ 2) - 1;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate(variable_index);

	const auto ref = 2 * X;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_13) {
	Polynomial X("x0");

	const auto p1 = (X ^ 2) + X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = 2 * X + 1;
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_14) {
	Polynomial X("x0");

	const auto p1 = (X ^ 2) + X + 1;

	constexpr size_t variable_index = 1;
	const auto result = p1.differentiate(variable_index);

	const auto ref = 0; 
	EXPECT_EQ(result, ref);
}
TEST(Polynomial, differentiate_15) {
	Polynomial X("x0");
	Polynomial Y("x1");

	const auto p1 = X * Y + X;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate(variable_index);

	const auto ref = Y + 1;
	EXPECT_EQ(result, ref);
}

//TEST(Polynomial, gradient_1) {
//	constexpr ushort space_dimension = 2;
//
//	const Polynomial x("x0");
//	const Polynomial y("x1");
//
//	const auto p = x + y + 1;
//	const auto result = p.gradient();
//
//	const Vector_Function<Polynomial,space_dimension> ref = { 1, 1 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Polynomial, gradient_2) {
//	constexpr ushort space_dimension = 2;
//
//	const Polynomial x("x0");
//	const Polynomial y("x1");
//
//	const auto p = (x ^ 2) + x * y + 1;
//	const auto result = p.gradient();
//
//	const Vector_Function<Polynomial, space_dimension> ref = { 2 * x + y, x };
//	EXPECT_EQ(result, ref);
//}
//TEST(Polynomial, gradient_3) {
//	const auto p = (X ^ 2) + X * Y + 1;
//	const auto result = p.gradient();
//
//	const Vector_Function<Polynomial, space_dimension> ref = { 2 * X + Y, X, 0 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Polynomial, gradient_4) {
//	const auto p = X * Y * Z;
//	const auto result = p.gradient();
//
//	const Vector_Function<Polynomial, space_dimension> ref = { Y * Z,  X * Z, X * Y };
//	EXPECT_EQ(result, ref);
//}
//
//
//
////Irrational Function
//TEST(IrrationalFunction, operator_call_1) {
//	constexpr size_t space_dimension = 1;
//
//	const auto x = Polynomial("x0");
//
//	const auto p = x;
//	const auto ir = p.root(0.5);
//
//	for (size_t i = 0; i < 10; ++i) {
//		const double val = 0.31 * i;
//		const std::vector<double> val_vec = { val };
//		const double result = ir(val_vec);
//		const double ref = std::pow(val, 0.5);
//		EXPECT_DOUBLE_EQ(result, ref);
//	}
//}
//TEST(IrrationalFunction, operator_call_2) {
//	constexpr size_t space_dimension = 1;
//
//	const auto x = Polynomial("x0");
//
//	const auto p = (x ^ 2) + x + 1;
//	const auto ir = p.root(0.5);
//
//	for (size_t i = 0; i < 10; ++i) {
//		const double val = 0.31 * i;
//		const std::vector<double> val_vec = { val };
//		const double result = ir(val_vec);
//		const double ref = std::pow(val * val + val + 1, 0.5);
//		EXPECT_DOUBLE_EQ(result, ref);
//	}
//}
//
//
//TEST(ms, is_natural_number_1) {
//	EXPECT_FALSE(ms::is_natural_number(3.1));
//}
//TEST(ms, is_natural_number_2) {
//	EXPECT_FALSE(ms::is_natural_number(-1));
//}









//TEST(Polynomial, OPERATOR_MULTIPLICATION_10) {
//	auto p1 = X + 1;
//	const auto p2 = p1^(0.5);
//	const auto result = p2 * p2;
//
//	auto ref = X + 1;
//	EXPECT_EQ(result, ref);
//}
//TEST(Polynomial, COMPLEX_OPERATION_2) {
//	auto p1 = (X ^ 2) + X + 1;
//	const auto p2 = p1 ^ 0.7;
//	const auto p3 = 5 * p2;
//
	//for (size_t i = 0; i < 10; ++i) {
	//	const double val = 0.31 * i;
	//	const MathVector val_vec = { val };
	//	const double result = p3(val_vec);
	//	const double ref = 5 * (std::pow(val * val + val + 1, 0.7));
	//	EXPECT_DOUBLE_EQ(result, ref);
//	}
//}
//TEST(Polynomial, COMPLEX_OPERATION_4) {
//	auto p1 = ms::sqrt(X + 1) * X + X + 1;
//	auto p2 = ms::sqrt(X + 2) * X;
//	const auto result = p1 * p2;
//
//	for (size_t i = 0; i < 10; ++i) {
//		const double val = 0.31 * i;
//		const MathVector v = { val };
//
//		const double ref = (std::sqrt(val + 1) * val + val + 1) * (std::sqrt(val + 2) * val);
//		EXPECT_DOUBLE_EQ(result(v), ref);
//	}
//}







//
//
//TEST(VECTORFUNCTION, OPERATOR_CALL1) {
//	auto p1 = (X ^ 2) + 3 * (X ^ 2) * Y + (Y ^ 3) + (Z ^ 2) - 6;
//	auto p2 = X + Y + Z - 3;
//	auto p3 = (Y ^ 2) * Z + X * Z - 2;
//	VectorFunction<auto> f = { p1,p2,p3 };
//	MathVector node = { 1,1,1 };
//	const auto result = f(node);
//
//	MathVector ref = { 0,0,0 };
//	EXPECT_EQ(result, ref);
//}
//
//
//TEST(VECTORFUNCTION, DIFFERENTIATE1) {
//	auto p1 = X + 1;
//	auto p2 = X + Y + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 1,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE2) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { Y,1,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE3) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { X,Z,1 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE4) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 1.5,0,0 };
//	EXPECT_EQ(f, ref);
//}
//TEST(VECTORFUNCTION, DIFFERENTIATE5) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate(variable_index);
//
//	VectorFunction<auto> ref = { 0.5,1,0 };
//	EXPECT_EQ(f, ref);
//}
//
//
//TEST(VECTORFUNCTION, CROSS_PRODUCT1) {
//	VectorFunction<auto> vf1 = { 1.5,0,0 };
//	VectorFunction<auto> vf2 = { 0.5,1,0 };
//	const auto result = vf1.cross_product(vf2);
//
//	VectorFunction<auto> ref = { 0,0,1.5 };
//	EXPECT_EQ(result, ref);
//}
//
//
//TEST(VECTORFUNCTION, L2_NORM1) {
//	VectorFunction<auto> vf = { 0,0,1.5 };
//	const auto result = vf.L2_norm();
//
//	auto ref = 1.5;
//	EXPECT_EQ(result, ref);
//}
