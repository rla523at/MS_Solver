#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Polynomial.h"

inline constexpr size_t space_dimension = 3;

#define X Polynomial<space_dimension>("x0")
#define Y Polynomial<space_dimension>("x1")
#define Z Polynomial<space_dimension>("x2")

GTEST_TEST(Polynomial, test_bed) {
}

GTEST_TEST(Polynomial, operator_addition_1) {
	const auto p1 = 3 * X + Y;
	const auto p2 = 5 * (X ^ 2);
	const auto p3 = Y ^ 2;
	const auto p4 = -7;
	const auto result = p1 + p2 + p3 + p4;

	const auto ref = 5 * (X ^ 2) + (Y ^ 2) + 3 * X + Y - 7;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_2) {
	const auto p1 = X + 1;
	const auto p2 = X + Y + 1;
	const auto result = p1 + p2;

	const auto ref = 2 * X + Y + 2;
	EXPECT_EQ(result, ref);	
}
GTEST_TEST(Polynomial, operator_addition_3) {
	const auto p1 = (X ^ 2) + X;
	const auto p2 = X + Y + 1;
	const auto result = p1 + p2;

	auto ref = (X ^ 2) + 2 * X + Y + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_4) {
	const Polynomial<space_dimension> p1 = 0.0;
	const Polynomial<space_dimension> p2 = 0.0;
	const auto result = p1 + p2;

	const Polynomial<space_dimension> ref = 0.0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_5) {
	const auto p1 = -12 * (X ^ 2) * Y + 45 * X;
	const auto p2 = 12 * (X ^ 2) * Y - 45 * X;
	const auto result = p1 + p2;

	const Polynomial<space_dimension> ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_6) {
	const auto p1 = (X ^ 2) + X;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = (X ^ 2) + 2 * X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_7) {
	const Polynomial<space_dimension> p1 = 0.0;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_8) {
	const Polynomial<space_dimension> p1 = 0.0;
	const auto p2 = X + 1;
	const auto result = p2 + p1;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_9) {
	const Polynomial<space_dimension> p1 = 3;
	const auto p2 = X + 1;
	const auto result = p1 + p2;

	const auto ref = X + 4;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_10) {
	const auto p1 = X + Y + Z + 5;
	const auto p2 = X + 1;
	const auto p3 = X + Y + 3;
	const auto p4 = X + Y + Z - 9;
	const auto result = p1 + p2 + p3 + p4;

	const auto ref = 4 * X + 3 * Y + 2 * Z;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_addition_11) {
	const auto p1 = X;
	const auto p2 = Y;
	const auto result = p1 + p2;

	const auto ref = X + Y;
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Polynomial, operator_substraction_1) {
	const auto p1 = X + Y;
	const auto p2 = X + 3 * Y;
	const auto result = p1 - p2;

	const auto ref = -2 * Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_2) {
	const auto p1 = X * Y;
	const auto p2 = X + 3 * Y;
	const auto result = p1 - p2;

	auto ref = X * Y - X - 3 * Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_3) {
	const auto p1 = X * Y;
	const auto p2 = X;
	const auto result = p1 - p2;

	const auto ref = X * Y - X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_4) {
	const Polynomial<space_dimension> p1 = 0;
	const auto p2 = X;
	const auto result = p1 - p2;

	const auto ref = -1 * X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_5) {
	const Polynomial<space_dimension> p1 = 0;
	const auto p2 = 3 * (X * Y) - 3 * ((X ^ 2) * Z);
	const auto result = p1 - p2;

	const auto ref = -3 * (X * Y) + 3 * ((X ^ 2) * Z);
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_6) {
	const auto p1 = (X - 1) * (Y - 1);
	const Polynomial<space_dimension> p2 = 1;
	const auto result = p1 - p2;

	const auto ref = X * Y - X - Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_7) {
	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = X + 1;
	const auto result = p1 - p2;

	const auto ref = X * Y - 2 * X - Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_substraction_8) {
	const auto p1 = (X - 1) * (Y - 1);
	const auto p2 = 2 * ((X - 1) ^ 2) + X + Y - 1;
	const auto result = p1 - p2;

	const auto ref = -2 * (X ^ 2) + X * Y + 2 * X - 2 * Y;
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Polynomial, operator_multiplication_1) {
	const auto p1 = X + 1;
	const auto p2 = X - 2;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) - X - 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_2) {
	const auto p1 = (X ^ 2) - 1;
	const auto p2 = (X ^ 2) + 2 * X + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) + 2 * (X ^ 3) - 2 * X - 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_3) {
	const auto p1 = (X ^ 2) + 2 * X + 1;
	const auto p2 = (X ^ 2) - 2 * X + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) - 2 * (X ^ 2) + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_4) {
	const Polynomial<space_dimension> p1 = 3;
	const Polynomial<space_dimension> p2 = 2;
	const auto result = p1 * p2;

	const Polynomial<space_dimension> ref = 6;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_5) {
	const Polynomial<space_dimension> p1 = 3;
	const auto p2 = X - 1;
	const auto result = p1 * p2;

	const auto ref = 3 * X - 3;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_6) {
	const auto p1 = X;
	const auto p2 = X + Y + 1;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + X * Y + X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_7) {
	const auto p1 = X + Y + 1;
	const auto p2 = X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + X * Y + X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_8) {
	const auto p1 = X + Y;
	const auto p2 = Y + X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) + 2 * X * Y + (Y ^ 2);
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_9) {
	const auto p1 = X + Y;
	const auto result = -1 * p1;

	const auto ref = -1 * X - 1 * Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_10) {
	const auto p1 = X + 1;
	const auto p2 = X + 1;
	const auto result = p1 * p2;

	auto ref = (X ^ 2) + 2 * X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_11) {
	const Polynomial<space_dimension> p1 = 0;
	const Polynomial<space_dimension> p2 = 0;
	const auto result = p1 * p2;

	const Polynomial<space_dimension> ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_12) {
	const auto p1 = 1;
	const auto p2 = X + 1;
	const auto result = p1 * p2;

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_13) {
	const auto p1 = X + 1;
	const Polynomial<space_dimension> p2 = 0;
	const auto result = p1 * p2;

	const Polynomial<space_dimension> ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_14) {
	const auto p1 = (X ^ 2) + X + 1;
	const auto result = 5 * p1;

	const auto ref = 5 * (X ^ 2) + 5 * X + 5;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_15) {
	const auto p1 = 3 * X;
	const auto p2 = 2 * Y;
	const auto result = p1 * p2;

	const auto ref = 6 * X * Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_16) {
	const auto p1 = 3 * X * Y;
	const auto p2 = (2 * X) ^ 2;
	const auto result = p1 * p2;

	const auto ref = 12 * (X ^ 3) * Y;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_17) {
	const auto p1 = X * (Y ^ 2) * (Z ^ 3);
	const auto p2 = X * (Y ^ 2);
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) * (Y ^ 4) * (Z ^ 3);
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_18) {
	const auto p1 = X * (Y ^ 4);
	const auto p2 = X * (Z ^ 3);
	const auto result = p1 * p2;

	const auto ref = (X ^ 2) * (Y ^ 4) * (Z ^ 3);
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_multiplication_19) {
	const auto p1 = X * (X + Y);
	const auto p2 = (Y + X) * X;
	const auto result = p1 * p2;

	const auto ref = (X ^ 4) + 2 * (X ^ 3) * Y + ((X * Y) ^ 2);
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Polynomial, operator_call_1) {
	const Polynomial<space_dimension> m = 1;
	const EuclideanVector values = { 1,2,3 };
	const auto result = m(values);

	constexpr double ref = 1.0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_call_2) {
	const Polynomial<space_dimension> m = 1;
	const EuclideanVector values = { 0,0,0 };
	const auto result = m(values);

	constexpr double ref = 1.0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, operator_call_3) {
	const auto p = 3 * X * (Y ^ 2) * (Z ^ 3);
	const EuclideanVector values = { 1,2,3 };
	const auto result = p(values);

	constexpr double ref2 = 324;
	EXPECT_EQ(result, ref2);
}
GTEST_TEST(Polynomial, operator_call_4) {
	const auto p = X * (Y ^ 2) * (Z ^ 3);
	const EuclideanVector values = { 1.84,2.789,3.487946 };
	const auto result = p(values);

	constexpr double ref3 = 6.073291260986822e+02;
	EXPECT_DOUBLE_EQ(result, ref3);
}

GTEST_TEST(Polynomial, operator_power_1) {
	const auto p = -0.125 * X + 0.125 * Y + 0.5;
	const auto result = p ^ 2;

	const auto ref = 0.015625 * (X ^ 2) - 0.03125 * X * Y + 0.015625 * (Y ^ 2) - 0.125 * X + 0.125 * Y + 0.25;
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Polynomial, complex_operation_1) {
	const auto p1 = X + 1;
	const auto p2 = X + 2;
	const auto p3 = X + 3;
	const auto p4 = p1 * p2 + p3;
	const auto p5 = p4 * p4;
	const auto result = 2 * (p5 + p1);

	auto ref = 2 * (X ^ 4) + 16 * (X ^ 3) + 52 * (X ^ 2) + 82 * X + 52;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, complex_operation_2) {
	auto p1 = X + 1;
	const auto p2 = p1 ^ 2;
	const auto p3 = p2 + X + 2;
	const auto p4 = p3 + 2 * X + 3;
	const auto result = 5 * p4;

	const auto ref = 5 * (X ^ 2) + 25 * X + 30;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, complex_operation_3) {
	auto p1 = X + Y;
	auto p2 = Y + X;
	auto p3 = X;
	const auto result = (p1 * p3) * (p2 * p3);

	auto ref = (X ^ 4) + 2 * (X ^ 3) * Y + ((X * Y) ^ 2);
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, complex_operation_4) {
	auto p1 = 0.25 * Y + 1.25;
	auto p2 = -0.25 * Y + -0.75;
	auto p3 = 0.25 * X + 0.25;
	auto p4 = -0.25 * X + 0.25;
	const auto result = (p1 * p4) - (p2 * p3);

	auto ref = -0.125 * X + 0.125 * Y + 0.5;
	EXPECT_EQ(result, ref);
}


GTEST_TEST(Polynomial, order_1) {
	const Polynomial<space_dimension> p = 0;
	const auto result = p.order();

	constexpr size_t ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, order_2) {
	const auto p = X;
	const auto result = p.order();

	constexpr size_t ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, order_3) {
	const auto p = (X * Y)^2;
	const auto result = p.order();

	constexpr size_t ref = 4;
	EXPECT_EQ(result, ref);
}

GTEST_TEST(Polynomial, differentiate_1) {
	const auto p = (X + 1) * (X - 1) + 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 2 * X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_2) {
	const auto p = (X - 1) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 2 * X - 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_3) {
	const auto p = (X - 1) ^ 2;

	constexpr size_t variable_index = 1;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 0;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_4) {
	const auto p1 = (X + 1) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = 2 * X + 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_5) {
	const auto p1 = (X + Y + 2) * (X + 1);

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = 2 * X + Y + 3;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_6) {
	const auto p1 = (X + Y + 2) * (X + 1);

	constexpr size_t variable_index = 1;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_7) {
	const auto p1 = (X + 1) * (Y + 1) + X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = Y + 2;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_8) {
	const auto p1 = ((X + 1) ^ 2) * (X + 2);

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = 3 * (X ^ 2) + 8 * X + 5;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_9) {
	const auto p = (2 * X + 3) ^ 2;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 8 * X + 12;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_10) {
	const auto p = ((2 * X + 3) ^ 2) * X;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 12 * (X ^ 2) + 24 * X + 9;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_11) {
	const auto  p = X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	Polynomial<space_dimension> ref = 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_12) {
	const auto p = (X ^ 2) - 1;

	constexpr size_t variable_index = 0;
	const auto result = p.differentiate<variable_index>();

	const auto ref = 2 * X;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_13) {
	const auto p1 = (X ^ 2) + X + 1;

	constexpr size_t variable_index = 0;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = 2 * X + 1;
	EXPECT_EQ(result, ref);
}
GTEST_TEST(Polynomial, differentiate_14) {
	const auto p1 = (X ^ 2) + X + 1;

	constexpr size_t variable_index = 1;
	const auto result = p1.differentiate<variable_index>();

	const auto ref = 0; 
	EXPECT_EQ(result, ref);
}

//GTEST_TEST(Polynomial, gradient_1) {
//	const auto p = X + Y + 1;
//	const auto result = p.gradient();
//
//	const VectorFunction<Polynomial> ref = { 1, 1 };
//	EXPECT_EQ(result, ref);
//}
//GTEST_TEST(Polynomial, gradient_2) {
//	const auto p = (X ^ 2) + X * Y + 1;
//	const auto result = p.gradient();
//
//	const VectorFunction<Polynomial> ref = { 2 * X + Y, X };
//	EXPECT_EQ(result, ref);
//}
//GTEST_TEST(Polynomial, gradient_3) {
//	const auto p = (X ^ 2) + X * Y + 1;
//	const size_t domain_dimension = 3;
//	const auto result = p.gradient(domain_dimension);
//
//	const VectorFunction<Polynomial> ref = { 2 * X + Y, X, 0 };
//	EXPECT_EQ(result, ref);
//}
//GTEST_TEST(Polynomial, gradient_4) {
//	const auto p = X * Y * Z;
//	const auto result = p.gradient();
//
//	const VectorFunction<Polynomial> ref = { Y * Z,  X * Z, X * Y };
//	EXPECT_EQ(result, ref);
//}
//
//
//
////Irrational Function
//GTEST_TEST(IrrationalFunction, operator_call_1) {
//	const auto p = X;
//	const auto ir = p.root(0.5);
//
//	for (size_t i = 0; i < 10; ++i) {
//		const double val = 0.31 * i;
//		const MathVector val_vec = { val };
//		const double result = ir(val_vec);
//		const double ref = std::pow(val, 0.5);
//		EXPECT_DOUBLE_EQ(result, ref);
//	}
//}
//GTEST_TEST(IrrationalFunction, operator_call_2) {
//	const auto p = (X ^ 2) + X + 1;
//	const auto ir = p.root(0.5);
//
//	for (size_t i = 0; i < 10; ++i) {
//		const double val = 0.31 * i;
//		const MathVector val_vec = { val };
//		const double result = ir(val_vec);
//		const double ref = std::pow(val * val + val + 1, 0.5);
//		EXPECT_DOUBLE_EQ(result, ref);
//	}
//}


//GTEST_TEST(Polynomial, OPERATOR_MULTIPLICATION_10) {
//	auto p1 = X + 1;
//	const auto p2 = p1^(0.5);
//	const auto result = p2 * p2;
//
//	auto ref = X + 1;
//	EXPECT_EQ(result, ref);
//}
//GTEST_TEST(Polynomial, COMPLEX_OPERATION_2) {
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
//GTEST_TEST(Polynomial, COMPLEX_OPERATION_4) {
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
//GTEST_TEST(VECTORFUNCTION, OPERATOR_CALL1) {
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
//GTEST_TEST(VECTORFUNCTION, DIFFERENTIATE1) {
//	auto p1 = X + 1;
//	auto p2 = X + Y + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate<variable_index>();
//
//	VectorFunction<auto> ref = { 1,1,0 };
//	EXPECT_EQ(f, ref);
//}
//GTEST_TEST(VECTORFUNCTION, DIFFERENTIATE2) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate<variable_index>();
//
//	VectorFunction<auto> ref = { Y,1,0 };
//	EXPECT_EQ(f, ref);
//}
//GTEST_TEST(VECTORFUNCTION, DIFFERENTIATE3) {
//	auto p1 = X * Y + 1;
//	auto p2 = X + Y * Z + Z - 3;
//	auto p3 = Y + Z;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate<variable_index>();
//
//	VectorFunction<auto> ref = { X,Z,1 };
//	EXPECT_EQ(f, ref);
//}
//GTEST_TEST(VECTORFUNCTION, DIFFERENTIATE4) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 0;
//	f.differentiate<variable_index>();
//
//	VectorFunction<auto> ref = { 1.5,0,0 };
//	EXPECT_EQ(f, ref);
//}
//GTEST_TEST(VECTORFUNCTION, DIFFERENTIATE5) {
//	auto p1 = 1.5 * X + 0.5 * Y + 3;
//	auto p2 = Y + 3;
//	auto p3 = 0;
//	VectorFunction<auto> f = { p1,p2,p3 };
//
//	constexpr size_t variable_index = 1;
//	f.differentiate<variable_index>();
//
//	VectorFunction<auto> ref = { 0.5,1,0 };
//	EXPECT_EQ(f, ref);
//}
//
//
//GTEST_TEST(VECTORFUNCTION, CROSS_PRODUCT1) {
//	VectorFunction<auto> vf1 = { 1.5,0,0 };
//	VectorFunction<auto> vf2 = { 0.5,1,0 };
//	const auto result = vf1.cross_product(vf2);
//
//	VectorFunction<auto> ref = { 0,0,1.5 };
//	EXPECT_EQ(result, ref);
//}
//
//
//GTEST_TEST(VECTORFUNCTION, L2_NORM1) {
//	VectorFunction<auto> vf = { 0,0,1.5 };
//	const auto result = vf.L2_norm();
//
//	auto ref = 1.5;
//	EXPECT_EQ(result, ref);
//}
