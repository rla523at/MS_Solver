#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Geometry.h"

TEST(Geometry, center_1) {
	const Figure fig = Figure::line;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const auto result = geometry.center_node();

	const Euclidean_Vector ref = { 1.5,1 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, center_2) {
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const auto result = geometry.center_node();

	const Euclidean_Vector ref = { 7.0 / 3.0, 4.0 / 3.0 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, center_3) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	
	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { 1,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	
	const auto result = geometry.center_node();

	const Euclidean_Vector ref = { 2,1.5 };
	EXPECT_EQ(result, ref);
}
