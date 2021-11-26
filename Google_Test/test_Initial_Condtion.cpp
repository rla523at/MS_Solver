#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Initial_Condition.h"

#include <random>

TEST(Constant1, calculate_solution_1) 
{
	Constant1 initial_condition;	
	Euclidean_Vector space_vector = { 1,1 };
	const auto result = initial_condition.calculate_solution(space_vector);

	const Euclidean_Vector ref = { 1.0 };
	EXPECT_EQ(result, ref);
}
TEST(Constant1, calculate_solution_2) 
{
	Constant1 initial_condition;
	Euclidean_Vector space_vector = { 1,1,1 };
	const auto result = initial_condition.calculate_solution(space_vector);

	const Euclidean_Vector ref = { 1.0 };
	EXPECT_EQ(result, ref);
}

TEST(Sine_Wave_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);

	Sine_Wave_2D initial_condition(1, 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		const Euclidean_Vector ref = { std::sin(2 * std::numbers::pi * space_vector[0]) * std::sin(2 * std::numbers::pi * space_vector[1]) };
		EXPECT_EQ(result, ref);
	}
}

TEST(Sine_Wave_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	
	Sine_Wave_3D initial_condition(1, 1, 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		const Euclidean_Vector ref = { std::sin(2 * std::numbers::pi * space_vector[0]) * std::sin(2 * std::numbers::pi * space_vector[1]) * std::sin(2 * std::numbers::pi * space_vector[2]) };
		EXPECT_EQ(result, ref);
	}
}

TEST(Square_Wave_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Square_Wave_2D initial_condition;

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector ref = { 0.0 };
		if (0.25 <= space_vector[0] && space_vector[0] <= 0.75 && 0.25 <= space_vector[1] && space_vector[1] <= 0.75)
		{
			ref = { 1.0 };
		}

		EXPECT_EQ(result, ref);
	}
}

TEST(Square_Wave_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Square_Wave_3D initial_condition;

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector ref = { 0.0 };
		if (0.25 <= space_vector[0] && space_vector[0] <= 0.75 && 0.25 <= space_vector[1] && space_vector[1] <= 0.75 && 0.25 <= space_vector[2] && space_vector[2] <= 0.75)
		{
			ref = { 1.0 };
		}

		EXPECT_EQ(result, ref);
	}
}

TEST(Circle_Wave_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Circle_Wave_2D initial_condition;

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector ref = { 0.0 };
		if ((space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) <= 0.25 * 0.25)
		{
			ref = { 1.0 };
		}

		EXPECT_EQ(result, ref);
	}
}

TEST(Circle_Wave_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Circle_Wave_3D initial_condition;

	for (ushort i = 0; i < num_iter; ++i)
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector ref = {0.0};
		if ((space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) + (space_vector[2] - 0.5) * (space_vector[2] - 0.5) <= 0.25 * 0.25)
		{
			ref = {1.0};
		}

		EXPECT_EQ(result, ref);
	}
}

TEST(Gaussian_Wave_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Gaussian_Wave_2D initial_condition;

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		const auto temp = (space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5);
		constexpr auto beta = 20.0;

		Euclidean_Vector ref = { std::exp(-beta * temp) };
		EXPECT_EQ(result, ref);
	}
}

TEST(Gaussian_Wave_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Gaussian_Wave_3D initial_condition;

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		const auto temp = (space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) + (space_vector[2] - 0.5) * (space_vector[2] - 0.5);
		constexpr auto beta = 20.0;
		Euclidean_Vector ref = { std::exp(-beta * temp) };

		EXPECT_EQ(result, ref);
	}
}

TEST(SOD_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	SOD_2D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0,0,c };
		Euclidean_Vector right = { 0.125,0,0,c * 0.1 };

		if (space_vector[0] <= 0.5)
		{
			EXPECT_EQ(result, left);
		}
		else
		{
			EXPECT_EQ(result, right);
		}
	}
}

TEST(SOD_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	SOD_3D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0,0,0,c };
		Euclidean_Vector right = { 0.125,0,0,0,c * 0.1 };

		if (space_vector[0] <= 0.5)
		{
			EXPECT_EQ(result, left);
		}
		else
		{
			EXPECT_EQ(result, right);
		}
	}
}

TEST(Modified_SOD_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Modified_SOD_2D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0.75,0,c + 0.5 * 0.75 * 0.75 };
		Euclidean_Vector right = { 0.125,0,0,c * 0.1 };

		if (space_vector[0] <= 0.3)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}

TEST(Modified_SOD_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Modified_SOD_3D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0.75,0,0,c + 0.5 * 0.75 * 0.75 };
		Euclidean_Vector right = { 0.125,0,0,0,c * 0.1 };

		if (space_vector[0] <= 0.3)
		{
			EXPECT_EQ(result, left);
		}
		else
		{
			EXPECT_EQ(result, right);
		}
	}
}

TEST(Shu_Osher_2D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Shu_Osher_2D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 3.857143, 3.857143 * 2.629369, 0, 10.333333 * c + 0.5 * 3.857143 * 2.629369 * 2.629369 };
		Euclidean_Vector right = { 1 + 0.2 * std::sin(16 * std::numbers::pi * space_vector[0]), 0, 0, c };

		if (space_vector[0] <= 0.125)
		{
			EXPECT_EQ(result, left);
		}
		else
		{
			EXPECT_EQ(result, right);
		}
	}
}

TEST(Shu_Osher_3D, calculate_solution_1) 
{
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	Shu_Osher_3D initial_condition;

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) 
	{
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = initial_condition.calculate_solution(space_vector);

		Euclidean_Vector left = { 3.857143, 3.857143 * 2.629369, 0,0, 10.333333 * c + 0.5 * 3.857143 * 2.629369 * 2.629369 };
		Euclidean_Vector right = { 1 + 0.2 * std::sin(16 * std::numbers::pi * space_vector[0]), 0, 0, 0, c };

		if (space_vector[0] <= 0.125)
		{
			EXPECT_EQ(result, left);
		}
		else
		{
			EXPECT_EQ(result, right);
		}
	}
}