#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Reconstruction_Method.h"

TEST(Reconstruction_Method, activation_function_1) {
	HardTanh f;

	Dynamic_Euclidean_Vector v = { -1,0,0.5,1,2,3 };
	v.apply(f);

	Dynamic_Euclidean_Vector ref = { 0,0,0.5,1,1,1 };
	EXPECT_EQ(v, ref);
}
TEST(Reconstruction_Method, activation_function_2) {
	ReLU f;

	Dynamic_Euclidean_Vector v = { -1,0,0.5,1,2,3 };
	v.apply(f);

	Dynamic_Euclidean_Vector ref = { 0,0,0.5,1,2,3 };
	EXPECT_EQ(v, ref);
}