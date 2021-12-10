#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/MBLAS.h"

#include <vector>

TEST(MSBLAS, x_plus_assign_y_1)
{
	std::vector<double> x = { 1,2,3 };
	std::vector<double> y = { 4,2,3 };

	ms::BLAS::x_plus_assign_y(static_cast<int>(x.size()), x.data(), y.data());

	std::vector<double> ref = { 5,4,6 };
	EXPECT_EQ(x, ref);
}