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
TEST(Geometry, volume_1) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 1,2 };
	const Euclidean_Vector n3 = { 2,2 };
	const Euclidean_Vector n4 = { 2,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 1;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_2) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { 1,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_3) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { -100,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 51;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_4) {
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.5;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_5) {
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1.524,1 };
	const Euclidean_Vector n2 = { 2,1.154 };
	const Euclidean_Vector n3 = { 4.47411,2 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.01084153;
	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
	EXPECT_NEAR(result, ref, 9.0E-16);
}
TEST(Geometry, volume_6) {
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.5, 1.257 };
	const Euclidean_Vector n3 = { 2.4874, 1.24 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.362569100000000;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Geometry, volume_7) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.0016, 1.257 };
	const Euclidean_Vector n3 = { 1.0017, 1.24 };
	const Euclidean_Vector n4 = { 1.001, 2.577 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };


	Geometry geometry(fig,fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 8.939999999999641e-04;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Geometry, volume_8) {
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.5, 1.257 };
	const Euclidean_Vector n3 = { 2.4874, 1.24 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };


	Geometry geometry(fig,fig_order, std::move(nodes));

	for (ushort integrand_order = 0; integrand_order < 16; ++integrand_order) {
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand_order);

		double result = 0.0;
		for (const auto weight : quadrature_rule.weights)
			result += weight;

		const auto ref = 0.362569100000000;
		//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
		EXPECT_NEAR(result, ref, 1.0E-15);
	}
}

TEST(Geometry, orthonormal_basis_1) {
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 3,1 };
	const Euclidean_Vector n3 = { 4,1 };
	const Euclidean_Vector n4 = { 1,3 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
		for (ushort j = 0; j <= i; ++j) {
			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);

			if (i == j)
				max_error = std::max(max_error, std::abs(1 - result));
			else
				max_error = std::max(max_error, std::abs(result));
		}
	}

	constexpr double allowable_error = 2.0E-9;
	EXPECT_LE(max_error, allowable_error);
}

TEST(Geometry, orthonormal_basis_2) {

	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 0,2 };
	const Euclidean_Vector n2 = { 2,0 };
	const Euclidean_Vector n3 = { 2,2 };
	const Euclidean_Vector n4 = { 0,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
		for (ushort j = 0; j <= i; ++j) {
			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);

			if (i == j)
				max_error = std::max(max_error, std::abs(1 - result));
			else
				max_error = std::max(max_error, std::abs(result));
		}
	}

	constexpr double allowable_error = 9.0E-11;
	EXPECT_LE(max_error, allowable_error);

}
TEST(Geometry, orthonormal_basis_3) {

	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 2.4874,1.257 };
	const Euclidean_Vector n3 = { 3.4874,1.24 };
	const Euclidean_Vector n4 = { 1,2.577 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
		for (ushort j = 0; j <= i; ++j) {
			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);

			if (i == j)
				max_error = std::max(max_error, std::abs(1 - result));
			else
				max_error = std::max(max_error, std::abs(result));
		}
	}

	constexpr double allowable_error = 9.0E-9;
	EXPECT_LE(max_error, allowable_error);
}
TEST(Geometry, orthonormal_basis_4) {

	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148 };
	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218 };
	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647 };
	const Euclidean_Vector n4 = { 0.3359655464243111, 0.2910404355350647 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
		for (ushort j = 0; j <= i; ++j) {
			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);

			if (i == j)
				max_error = std::max(max_error, std::abs(1 - result));
			else
				max_error = std::max(max_error, std::abs(result));
		}
	}

	constexpr double allowable_error = 9.0E-13;
	EXPECT_LE(max_error, allowable_error);

}
TEST(Geometry, orthonormal_basis_5) {

	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148 };
	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218 };
	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
		for (ushort j = 0; j <= i; ++j) {
			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);

			if (i == j)
				max_error = std::max(max_error, std::abs(1 - result));
			else
				max_error = std::max(max_error, std::abs(result));
		}
	}

	constexpr double allowable_error = 9.0E-13;
	EXPECT_LE(max_error, allowable_error);
}

