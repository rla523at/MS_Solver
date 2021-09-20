#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Initial_Condition.h"

#include <random>


TEST(Constant1, calculate_solution_1) {
	constexpr ushort space_dimension = 2;

	Euclidean_Vector space_vector = { 1,1 };
	const auto result = Constant1<space_dimension>::calculate_solution(space_vector);

	const Euclidean_Vector ref = 1.0;
	EXPECT_EQ(result, ref);
}
TEST(Constant1, calculate_solution_2) {
	constexpr ushort space_dimension = 3;

	Euclidean_Vector space_vector = { 1,1,1 };
	const auto result = Constant1<space_dimension>::calculate_solution(space_vector);

	const Euclidean_Vector ref = 1.0;
	EXPECT_EQ(result, ref);
}

TEST(Sine_Wave, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 5;

	Sine_Wave<space_dimension>::initialize({ 1,1 });
	
	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Sine_Wave<space_dimension>::calculate_solution(space_vector);

		const Euclidean_Vector ref = std::sin(2 * std::numbers::pi * space_vector[0]) * std::sin(2 * std::numbers::pi * space_vector[1]);
		EXPECT_EQ(result, ref);
	}
}
TEST(Sine_Wave, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 5;

	Sine_Wave<space_dimension>::initialize({ 1,1,1 });

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-99, 99);
	
	for (ushort i = 0; i < num_iter; ++i) {

	Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
	const auto result = Sine_Wave<space_dimension>::calculate_solution(space_vector);

	const Euclidean_Vector ref = std::sin(2 * std::numbers::pi * space_vector[0]) * std::sin(2 * std::numbers::pi * space_vector[1]) * std::sin(2 * std::numbers::pi * space_vector[2]);
	EXPECT_EQ(result, ref);
	}
}

TEST(Square_Wave, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Square_Wave<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector ref = 0.0;
		if (0.25 <= space_vector[0] && space_vector[0] <= 0.75 && 0.25 <= space_vector[1] && space_vector[1] <= 0.75)
			ref = 1.0;

		EXPECT_EQ(result, ref);
	}
}
TEST(Square_Wave, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {

		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = Square_Wave<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector ref = 0.0;
		if (0.25 <= space_vector[0] && space_vector[0] <= 0.75 && 0.25 <= space_vector[1] && space_vector[1] <= 0.75 && 0.25 <= space_vector[2] && space_vector[2] <= 0.75)
			ref = 1.0;

		EXPECT_EQ(result, ref);
	}
}

TEST(Circle_Wave, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Circle_Wave<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector ref = 0.0;
		if ((space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) <= 0.25 * 0.25)
			ref = 1.0;

		EXPECT_EQ(result, ref);
	}
}
TEST(Circle_Wave, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {

		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = Circle_Wave<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector ref = 0.0;
		if ((space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) + (space_vector[2] - 0.5) * (space_vector[2] - 0.5) <= 0.25 * 0.25)
			ref = 1.0;
		
		EXPECT_EQ(result, ref);
	}
}

TEST(Gaussian_Wave, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Gaussian_Wave<space_dimension>::calculate_solution(space_vector);
		
		const auto temp = (space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5);
		constexpr auto beta = 20.0;
				
		Euclidean_Vector ref = std::exp(-beta * temp);
		EXPECT_EQ(result, ref);
	}
}
TEST(Gaussian_Wave, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 5;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	for (ushort i = 0; i < num_iter; ++i) {

		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = Gaussian_Wave<space_dimension>::calculate_solution(space_vector);

		const auto temp = (space_vector[0] - 0.5) * (space_vector[0] - 0.5) + (space_vector[1] - 0.5) * (space_vector[1] - 0.5) + (space_vector[2] - 0.5) * (space_vector[2] - 0.5);
		constexpr auto beta = 20.0;
		Euclidean_Vector ref = std::exp(-beta * temp);

		EXPECT_EQ(result, ref);
	}
}

TEST(SOD, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = SOD<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left  = { 1.0,0,0,c };
		Euclidean_Vector right = { 0.125,0,0,c * 0.1 };

		if (space_vector[0] < 0.5)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}
TEST(SOD, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = SOD<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0,0,0,c };
		Euclidean_Vector right = { 0.125,0,0,0,c * 0.1 };

		if (space_vector[0] < 0.5)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}

TEST(Modified_SOD, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Modified_SOD<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0.75,0,c + 0.5 * 0.75 * 0.75 };
		Euclidean_Vector right = { 0.125,0,0,c * 0.1 };

		if (space_vector[0] < 0.3)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}
TEST(Modified_SOD, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(0, 1);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = Modified_SOD<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left = { 1.0,0.75,0,0,c + 0.5 * 0.75 * 0.75 };
		Euclidean_Vector right = { 0.125,0,0,0,c * 0.1 };

		if (space_vector[0] < 0.3)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}

TEST(Shu_Osher, calculate_solution_1) {
	constexpr ushort space_dimension = 2;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-5, 5);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen) };
		const auto result = Shu_Osher<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left = { 3.857143, 3.857143 * 2.629369, 0, 10.333333 * c + 0.5 * 3.857143 * 2.629369 * 2.629369 };
		Euclidean_Vector right = { 1 + 0.2 * std::sin(5 * space_vector[0]), 0, 0, c };

		if (space_vector[0] < -4)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}
TEST(Shu_Osher, calculate_solution_2) {
	constexpr ushort space_dimension = 3;
	constexpr ushort num_iter = 10;

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis1(-5, 5);

	constexpr auto gamma = 1.4;
	constexpr auto c = 1 / (gamma - 1);

	for (ushort i = 0; i < num_iter; ++i) {
		Euclidean_Vector space_vector = { dis1(gen),dis1(gen),dis1(gen) };
		const auto result = Shu_Osher<space_dimension>::calculate_solution(space_vector);

		Euclidean_Vector left = { 3.857143, 3.857143 * 2.629369, 0,0, 10.333333 * c + 0.5 * 3.857143 * 2.629369 * 2.629369 };
		Euclidean_Vector right = { 1 + 0.2 * std::sin(5 * space_vector[0]), 0, 0, 0, c };

		if (space_vector[0] < -4)
			EXPECT_EQ(result, left);
		else
			EXPECT_EQ(result, right);
	}
}