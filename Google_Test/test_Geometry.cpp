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

	const auto result = geometry.center_point();

	const Euclidean_Vector ref = { 1.5,1 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, center_2) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const auto result = geometry.center_point();

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
	
	const auto result = geometry.center_point();

	const Euclidean_Vector ref = { 2,1.5 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, center_4) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1,2 };
	const Euclidean_Vector n2 = { 2,1,3 };
	const Euclidean_Vector n3 = { 4,2,4 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const auto result = geometry.center_point();

	const Euclidean_Vector ref = { 7.0 / 3.0, 4.0 / 3.0, 3.0 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, change_points_1)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const Euclidean_Vector n4 = { 5,4 };
	const Euclidean_Vector n5 = { 6,3 };
	const Euclidean_Vector n6 = { 7,2 };
	std::vector<Euclidean_Vector> nodes2 = { n4,n5,n6 };

	geometry.change_points(std::move(nodes2));

	const auto result = geometry.center_point();

	const Euclidean_Vector ref = { 6, 3 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, change_points_2)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const Euclidean_Vector n4 = { 1.524,1 };
	const Euclidean_Vector n5 = { 2,1.154 };
	const Euclidean_Vector n6 = { 4.47411,2 };
	std::vector<Euclidean_Vector> nodes2 = { n4,n5,n6 };

	geometry.change_points(std::move(nodes2));

	const auto result = geometry.volume();

	const auto ref = 0.01084153;
	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
	EXPECT_NEAR(result, ref, 9.0E-16);
}
TEST(Geometry, face_geometries_1)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto result = geometry.face_geometries();

	const Figure f_fig = Figure::line;

	const Euclidean_Vector f1_n1 = { 1,1 };
	const Euclidean_Vector f1_n2 = { 2,1 };
	std::vector<Euclidean_Vector> f1_nodes = { f1_n1,f1_n2 };
	Geometry f1_geometry(f_fig, fig_order, std::move(f1_nodes));

	const Euclidean_Vector f2_n1 = { 2,1 };
	const Euclidean_Vector f2_n2 = { 4,2 };
	std::vector<Euclidean_Vector> f2_nodes = { f2_n1,f2_n2 };
	Geometry f2_geometry(f_fig, fig_order, std::move(f2_nodes));

	const Euclidean_Vector f3_n1 = { 4,2 };
	const Euclidean_Vector f3_n2 = { 1,1 };
	std::vector<Euclidean_Vector> f3_nodes = { f3_n1,f3_n2 };
	Geometry f3_geometry(f_fig, fig_order, std::move(f3_nodes));

	//const std::vector<Geometry> ref = { std::move(f1_geometry),std::move(f2_geometry),std::move(f3_geometry) };
	EXPECT_EQ(result[0], f1_geometry);
	EXPECT_EQ(result[1], f2_geometry);
	EXPECT_EQ(result[2], f3_geometry);
}
TEST(Geometry, normalized_normal_vector_1)
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };

	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto normal = geometry.normalized_normal_vector(geometry.center_point());

	const Euclidean_Vector ref_direct = { 0,1 };

	const auto tangent = n2 - n1;
	const auto normality = tangent.inner_product(normal);
	const double size = normal.L2_norm();
	const double direction = normal.inner_product(ref_direct);

	const double normality_ref = 0;
	const double size_ref = 1;
	const double direction_ref = 1;

	EXPECT_EQ(normality, normality_ref);
	EXPECT_EQ(size, size_ref);
	EXPECT_EQ(direction, direction_ref);
}
TEST(Geometry, normalized_normal_vector_2)
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };


	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto normal = geometry.normalized_normal_vector(geometry.center_point());

	const auto tangent = n2 - n1;
	const auto normality = tangent.inner_product(normal);
	const double size = normal.L2_norm();

	const double normality_ref = 0;
	const double size_ref = 1;

	EXPECT_EQ(normality, normality_ref);
	EXPECT_DOUBLE_EQ(size, size_ref);
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
	for (ushort i = 0; i < orthonormal_basis.size(); ++i) {
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
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 0,2 };
	const Euclidean_Vector n2 = { 2,0 };
	const Euclidean_Vector n3 = { 2,2 };
	const Euclidean_Vector n4 = { 0,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.size(); ++i) {
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
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 2.4874,1.257 };
	const Euclidean_Vector n3 = { 3.4874,1.24 };
	const Euclidean_Vector n4 = { 1,2.577 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.size(); ++i) {
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
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148 };
	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218 };
	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647 };
	const Euclidean_Vector n4 = { 0.3359655464243111, 0.2910404355350647 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.size(); ++i) {
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
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148 };
	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218 };
	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	constexpr auto polynomial_order = 5;
	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function(polynomial_order);

	double max_error = 0.0;
	for (ushort i = 0; i < orthonormal_basis.size(); ++i) {
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
TEST(Geometry, projected_volume_1)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1.524,1 };
	const Euclidean_Vector n2 = { 2,1.154 };
	const Euclidean_Vector n3 = { 4.47411,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };


	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto result = geometry.projected_volumes();

	const std::vector<double> ref = { 1, 4.47411 - 1.524 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, projected_volume_2)
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { 1,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto result = geometry.projected_volumes();

	const std::vector<double> ref = { 1, 3 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, projected_volume_3)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1,2 };
	const Euclidean_Vector n2 = { 2,1,2 };
	const Euclidean_Vector n3 = { 4,2,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));

	const auto result = geometry.projected_volumes();

	const std::vector<double> ref = { 0.0, 0.0, 0.5 };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, set_of_face_points1)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1.524,1 };
	const Euclidean_Vector n2 = { 2,1.154 };
	const Euclidean_Vector n3 = { 4.47411,2 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };


	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto result = geometry.set_of_face_points();

	const std::vector<std::vector<Euclidean_Vector>> ref = { {n1,n2},{n2,n3},{n3,n1} };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, set_of_face_points2)
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	Euclidean_Vector n1 = { 1,1 };
	Euclidean_Vector n2 = { 2,1 };
	Euclidean_Vector n3 = { 4,2 };
	Euclidean_Vector n4 = { 1,2 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(std::move(ref_geo), std::move(nodes));
	const auto result = geometry.set_of_face_points();

	const std::vector<std::vector<Euclidean_Vector>> ref = { {n1,n2},{n2,n3},{n3,n4}, {n4,n1} };
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_1) 
{
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
TEST(Geometry, volume_2) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { 1,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 2;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_3) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { -100,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 51;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_4)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.5;
	EXPECT_EQ(result, ref);
}
TEST(Geometry, volume_5)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1.524,1 };
	const Euclidean_Vector n2 = { 2,1.154 };
	const Euclidean_Vector n3 = { 4.47411,2 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.01084153;
	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
	EXPECT_NEAR(result, ref, 9.0E-16);
}
TEST(Geometry, volume_6)
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.5, 1.257 };
	const Euclidean_Vector n3 = { 2.4874, 1.24 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };

	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 0.362569100000000;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Geometry, volume_7)
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.0016, 1.257 };
	const Euclidean_Vector n3 = { 1.0017, 1.24 };
	const Euclidean_Vector n4 = { 1.001, 2.577 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };


	Geometry geometry(fig, fig_order, std::move(nodes));
	const auto result = geometry.volume();

	const auto ref = 8.939999999999641e-04;
	EXPECT_DOUBLE_EQ(result, ref);
}
TEST(Geometry, volume_8) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	const Euclidean_Vector n1 = { 1,2 };
	const Euclidean_Vector n2 = { 1.5, 1.257 };
	const Euclidean_Vector n3 = { 2.4874, 1.24 };

	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };


	Geometry geometry(fig, fig_order, std::move(nodes));

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
//TEST(Geometry, sub_simplex1) {
//
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(std::move(ref_geo), std::move(nodes));
//
//	const auto result = geometry.sub_simplexgeometries();
//
//	const auto ref_fig = Figure::triangle;
//	const ReferenceGeometry simplexref_geometry(ref_fig, fig_order);
//
//	std::vector<Euclidean_Vector> simplexnodes1 = { n1,n2,n4 };
//	std::vector<Euclidean_Vector> simplexnodes2 = { n2,n3,n1 };
//	std::vector<Euclidean_Vector> simplexnodes3 = { n3,n4,n2 };
//	std::vector<Euclidean_Vector> simplexnodes4 = { n4,n1,n3 };
//
//	const Geometry simplex1(simplexref_geometry, std::move(simplexnodes1));
//	const Geometry simplex2(simplexref_geometry, std::move(simplexnodes2));
//	const Geometry simplex3(simplexref_geometry, std::move(simplexnodes3));
//	const Geometry simplex4(simplexref_geometry, std::move(simplexnodes4));
//
//	std::vector<Geometry> ref = { simplex1,simplex2,simplex3,simplex4 };
//	EXPECT_EQ(result, ref);
//}


//TEST(Geometry, projected_volume_3) 
//{
//	constexpr ushort space_dimension = 3;
//
//	const auto fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1,0 };
//	const Euclidean_Vector n2 = { 2,1,0 };
//	const Euclidean_Vector n3 = { 2,2,0 };
//	const Euclidean_Vector n4 = { 1,2,0 };
//	const Euclidean_Vector n5 = { 1,1,2 };
//	const Euclidean_Vector n6 = { 2,1,2 };
//	const Euclidean_Vector n7 = { 2,2,2 };
//	const Euclidean_Vector n8 = { 1,2,2 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4,n5,n6,n7,n8 };
//
//	Geometry geometry(std::move(ref_geo), std::move(nodes));
//	const auto result = geometry.projected_volume();
//
//	const std::array<double, space_dimension> ref = { 2,2,1 };
//	EXPECT_EQ(result, ref);
//}




//TEST(Geometry, is_axis_parallel_1) 
//{
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0,0 };
//	const Euclidean_Vector n2 = { 0,1 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2 };
//	Geometry geometry(std::move(ref_geo), std::move(nodes));
//
//	const Euclidean_Vector n3 = { 1,1 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector> nodes2 = { n3,n4 };
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	EXPECT_FALSE(geometry.can_be_periodic_pair(geometry2));
//}
//TEST(Geometry, is_on_same_axis_1) {
//	constexpr ushort space_dimension = 3;
//
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148, 0 };
//	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218, 0 };
//	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647, 0 };
//	std::vector<Euclidean_Vector> nodes1 = { n1,n2,n3 };
//
//	const Euclidean_Vector n4 = { 0.3635520579711813, 0.1973431147402148, 0 };
//	const Euclidean_Vector n5 = { 0.3512301560533574, 0.2184608229801218, 0 };
//	const Euclidean_Vector n6 = { 0.3309655464243111, 0.2010404355350647, 0 };
//	std::vector<Euclidean_Vector> nodes2 = { n4,n5,n6 };
//
//	Geometry geometry1(ref_geometry, std::move(nodes1));
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	EXPECT_TRUE(geometry1.is_on_same_axis(geometry2));
//}
//TEST(Geometry, is_on_same_axis_2) {
//	constexpr ushort space_dimension = 3;
//
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148, 1 };
//	const Euclidean_Vector n2 = { 0.3512301560533574, 0.2973431147402148, 2 };
//	const Euclidean_Vector n3 = { 0.3309655464243111, 0.2973431147402148, 3 };
//	std::vector<Euclidean_Vector> nodes1 = { n1,n2,n3 };
//
//	const Euclidean_Vector n4 = { 0.3635520579711813, 0.2973431147402148, 3 };
//	const Euclidean_Vector n5 = { 0.3512301560533574, 0.2973431147402148, 2 };
//	const Euclidean_Vector n6 = { 0.3309655464243111, 0.2973431147402148, 1 };
//	std::vector<Euclidean_Vector> nodes2 = { n4,n5,n6 };
//
//	Geometry geometry1(ref_geometry, std::move(nodes1));
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	EXPECT_TRUE(geometry1.is_on_same_axis(geometry2));
//}