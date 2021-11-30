#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Configuration.h"

TEST(Configuration, get1)
{
	const auto configuration_file_path = "RSC/configuration.dat";
	Configuration configuration(configuration_file_path);

	const auto result = configuration.get("grid_file_path");
	
	const auto ref = "RSC/Grid/2D/RQ10.msh";
	EXPECT_EQ(result, ref);
}
TEST(Configuration, get2)
{
	const auto configuration_file_path = "RSC/configuration.dat";
	Configuration configuration(configuration_file_path);

	const auto result = configuration.get<short>("space_dimension");

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}