#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Element.h"


//TEST(ReferenceGeometry, scale_function_1) {
///
//	const Figure fig = Figure::quadrilateral;
//	const order fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2 };
//	Euclidean_Vector p2 = { 3,1 };
//	Euclidean_Vector p3 = { 4,1 };
//	Euclidean_Vector p4 = { 1,3 };
//	std::vector<Euclidean_Vector> pv = { p1,p2,p3,p4 };
//
//	const auto mapping_function = ref_geometry.mapping_function(pv);
//	const auto result = ref_geometry.scale_function(mapping_function);
//
//	const auto x = Polynomial("x0");
//	const auto y = Polynomial("x1");
//	Irrational_Function ref = (0.25 * y + 1.25) * (-0.25 * x + 0.25) - (0.25 * x + 0.25) * (-0.25 * y - 0.75);
//	EXPECT_EQ(result, ref);
//}

//TEST(Geometry, volume_1) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 1,2 };
//	const Euclidean_Vector n3 = { 2,2 };
//	const Euclidean_Vector n4 = { 2,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 1;
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, volume_2) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 2;
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, volume_3) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { -100,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 51;
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, volume_4) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 0.5;
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, volume_5) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1.524,1 };
//	const Euclidean_Vector n2 = { 2,1.154 };
//	const Euclidean_Vector n3 = { 4.47411,2 };
//
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 0.01084153;
//	//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
//	EXPECT_NEAR(result, ref, 9.0E-16);
//}
//TEST(Geometry, volume_6) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,2 };	
//	const Euclidean_Vector n2 = { 1.5, 1.257 };
//	const Euclidean_Vector n3 = { 2.4874, 1.24 };
//
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 0.362569100000000;
//	EXPECT_DOUBLE_EQ(result, ref);
//}
//TEST(Geometry, volume_7) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,2 };
//	const Euclidean_Vector n2 = { 1.0016, 1.257 };
//	const Euclidean_Vector n3 = { 1.0017, 1.24 };
//	const Euclidean_Vector n4 = { 1.001, 2.577};
//
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.volume();
//
//	const auto ref = 8.939999999999641e-04;
//	EXPECT_DOUBLE_EQ(result, ref);
//}
//
//
//TEST(Geometry, quadrature_weight_1) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,2 };
//	const Euclidean_Vector n2 = { 1.5, 1.257 };
//	const Euclidean_Vector n3 = { 2.4874, 1.24 };
//
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	
//	for (ushort integrand_order = 0; integrand_order < 16; ++integrand_order) {
//		const auto quadrature_rule = geometry.get_quadrature_rule(integrand_order);
//
//		double result = 0.0;
//		for (const auto weight : quadrature_rule.weights)
//			result += weight;
//
//		const auto ref = 0.362569100000000;
//		//EXPECT_DOUBLE_EQ(result, ref); //suffer by round off error
//		EXPECT_NEAR(result, ref, 1.0E-15);
//	}
//}
//
//
//TEST(Geometry, faces_nodes1) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1.524,1 };
//	const Euclidean_Vector n2 = { 2,1.154 };
//	const Euclidean_Vector n3 = { 4.47411,2 };
//
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.set_of_face_nodes();
//
//	const std::vector<std::vector<Euclidean_Vector>> ref = { {n1,n2},{n2,n3},{n3,n1} };
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, faces_nodes2) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//
//	Euclidean_Vector n1 = { 1,1 };
//	Euclidean_Vector n2 = { 2,1 };
//	Euclidean_Vector n3 = { 4,2 };
//	Euclidean_Vector n4 = { 1,2 };
//
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.set_of_face_nodes();
//
//	const std::vector<std::vector<Euclidean_Vector>> ref = { {n1,n2},{n2,n3},{n3,n4}, {n4,n1} };
//	EXPECT_EQ(result, ref);
//}
//
//
//TEST(Geometry, sub_simplex1) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	
//	const auto result = geometry.sub_simplexgeometries();
//
//	const auto ref_fig = Figure::triangle;
//	const ReferenceGeometry simplexref_geometry(ref_fig, fig_order);
//
//	std::vector<Euclidean_Vector<2>> simplexnodes1 = { n1,n2,n4 };
//	std::vector<Euclidean_Vector<2>> simplexnodes2 = { n2,n3,n1 };
//	std::vector<Euclidean_Vector<2>> simplexnodes3 = { n3,n4,n2 };
//	std::vector<Euclidean_Vector<2>> simplexnodes4 = { n4,n1,n3 };
//
//	const Geometry simplex1(simplexref_geometry, std::move(simplexnodes1));
//	const Geometry simplex2(simplexref_geometry, std::move(simplexnodes2));
//	const Geometry simplex3(simplexref_geometry, std::move(simplexnodes3));
//	const Geometry simplex4(simplexref_geometry, std::move(simplexnodes4));
//
//	std::vector<Geometry> ref = { simplex1,simplex2,simplex3,simplex4 };
//	EXPECT_EQ(result, ref);
//}
//
//TEST(Geometry, projected_volume_1) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1.524,1 };
//	const Euclidean_Vector n2 = { 2,1.154 };
//	const Euclidean_Vector n3 = { 4.47411,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.projected_volume();
//
//	const std::array<double, 2> ref = { 1, 4.47411 - 1.524 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, projected_volume_2) {
///
//
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.projected_volume();
//
//	const std::array<double, 2> ref = { 1, 3 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Geometry, projected_volume_3) {
//	constexpr ushort space_dimension = 3;
//
//	const auto fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
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
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.projected_volume();
//
//	const std::array<double, space_dimension> ref = { 2,2,1 };
//	EXPECT_EQ(result, ref);
//}
//

//TEST(Geometry, normalized_normal_vector_1) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto normal = geometry.normalized_normal_vector(geometry.center_node());
//
//	const auto tangent = n2 - n1;
//	const auto normality = tangent.inner_product(normal);
//	const double size = normal.L2_norm();
//	const double direction = normal.inner_product({ 0,1 });
//
//
//	const double normality_ref = 0;
//	const double size_ref = 1;
//	const double direction_ref = 1;
//
//	EXPECT_EQ(normality, normality_ref);
//	EXPECT_EQ(size, size_ref);
//	EXPECT_EQ(direction, direction_ref);
//}
//TEST(Geometry, normalized_normal_vector_2) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto normal = geometry.normalized_normal_vector(geometry.center_node());
//
//	const auto tangent = n2 - n1;
//	const auto normality = tangent.inner_product(normal);
//	const double size = normal.L2_norm();
//
//
//	const double normality_ref = 0;
//	const double size_ref = 1;
//
//	EXPECT_EQ(normality, normality_ref);
//	EXPECT_DOUBLE_EQ(size, size_ref);
//}
//
//TEST(Geometry, is_axis_parallel_1) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0,0 };
//	const Euclidean_Vector n2 = { 0,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	const Euclidean_Vector n3 = { 1,1 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes2 = { n3,n4 };
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	EXPECT_FALSE(geometry.can_be_periodic_pair(geometry2));
//}
//
//TEST(Geometry, orthonormal_basis_1) {
//	
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,2 };
//	const Euclidean_Vector n2 = { 3,1 };
//	const Euclidean_Vector n3 = { 4,1 };
//	const Euclidean_Vector n4 = { 1,3 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	constexpr auto polynomial_order = 5;
//	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function<polynomial_order>();
//
//	double max_error = 0.0;
//	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
//		for (ushort j = 0; j <= i; ++j) {
//			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);
//
//			if (i == j)
//				max_error = std::max(max_error, std::abs(1 - result));
//			else
//				max_error = std::max(max_error, std::abs(result));
//		}
//	}
//
//	constexpr double allowable_error = 2.0E-9;
//	EXPECT_LE(max_error, allowable_error);
//}
//
//TEST(Geometry, orthonormal_basis_2) {
//	
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0,2 };
//	const Euclidean_Vector n2 = { 2,0 };
//	const Euclidean_Vector n3 = { 2,2 };
//	const Euclidean_Vector n4 = { 0,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	constexpr auto polynomial_order = 5;
//	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function<polynomial_order>();
//
//	double max_error = 0.0;
//	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
//		for (ushort j = 0; j <= i; ++j) {
//			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);
//
//			if (i == j)
//				max_error = std::max(max_error, std::abs(1 - result));
//			else
//				max_error = std::max(max_error, std::abs(result));
//		}
//	}
//
//	constexpr double allowable_error = 9.0E-11;
//	EXPECT_LE(max_error, allowable_error);
//	
//}
//TEST(Geometry, orthonormal_basis_3) {
//	
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,2 };
//	const Euclidean_Vector n2 = { 2.4874,1.257 };
//	const Euclidean_Vector n3 = { 3.4874,1.24 };
//	const Euclidean_Vector n4 = { 1,2.577 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	constexpr auto polynomial_order = 5;
//	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function<polynomial_order>();
//
//	double max_error = 0.0;
//	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
//		for (ushort j = 0; j <= i; ++j) {
//			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);
//
//			if (i == j)
//				max_error = std::max(max_error, std::abs(1 - result));
//			else
//				max_error = std::max(max_error, std::abs(result));
//		}
//	}
//
//	constexpr double allowable_error = 9.0E-9;
//	EXPECT_LE(max_error, allowable_error);
//}
//TEST(Geometry, orthonormal_basis_4) {
//	
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148};
//	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218};
//	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647};
//	const Euclidean_Vector n4 = { 0.3359655464243111, 0.2910404355350647};
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	constexpr auto polynomial_order = 5;
//	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function<polynomial_order>();
//
//	double max_error = 0.0;
//	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
//		for (ushort j = 0; j <= i; ++j) {
//			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);
//
//			if (i == j)
//				max_error = std::max(max_error, std::abs(1 - result));
//			else
//				max_error = std::max(max_error, std::abs(result));
//		}
//	}
//
//	constexpr double allowable_error = 9.0E-13;
//	EXPECT_LE(max_error, allowable_error);
//	
//}
//TEST(Geometry, orthonormal_basis_5) {
//	
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 0.3635520579711813, 0.2973431147402148 };
//	const Euclidean_Vector n2 = { 0.3512301560533574, 0.3184608229801218 };
//	const Euclidean_Vector n3 = { 0.3309655464243111, 0.3010404355350647 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	constexpr auto polynomial_order = 5;
//	const auto orthonormal_basis = geometry.orthonormal_basis_vector_function<polynomial_order>();
//
//	double max_error = 0.0;
//	for (ushort i = 0; i < orthonormal_basis.range_dimension(); ++i) {
//		for (ushort j = 0; j <= i; ++j) {
//			const auto result = ms::inner_product(orthonormal_basis[i], orthonormal_basis[j], geometry);
//
//			if (i == j)
//				max_error = std::max(max_error, std::abs(1 - result));
//			else
//				max_error = std::max(max_error, std::abs(result));
//		}
//	}
//
//	constexpr double allowable_error = 9.0E-13;
//	EXPECT_LE(max_error, allowable_error);
//}
//
//TEST(Geometry, is_on_same_axis_1) {
//	constexpr ushort space_dimension = 3;
//
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
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
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
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
//
//TEST(Element, vertex_node_indexes_1) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	ElementType element_type = ElementType::cell;
//	std::vector<uint> indexes = { 5,6,7,8 };
//
//	Element element(element_type, std::move(geometry), std::move(indexes));
//	const auto result = element.vertex_node_indexes();
//
//	const std::vector<uint> ref = { 5,6,7,8 };
//	EXPECT_EQ(ref, result);
//}
//
//
//TEST(Element, set_of_face_node_indexes) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	ElementType element_type = ElementType::cell;
//	std::vector<uint> indexes = { 5,6,7,8 };
//
//	Element element(element_type, std::move(geometry), std::move(indexes));
//	const auto result = element.set_of_face_node_indexes();
//
//	const std::vector<std::vector<uint>> ref = { {5,6},{6,7},{7,8},{8,5} };
//	EXPECT_EQ(ref, result);
//}
//
//
//TEST(Element, set_of_face_vertex_node_indexes) {
///
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	ElementType element_type = ElementType::cell;
//	std::vector<uint> indexes = { 5,6,7,8 };
//
//	Element element(element_type, std::move(geometry), std::move(indexes));
//	const auto result = element.set_of_face_vertex_node_indexes();
//
//	const std::vector<std::vector<uint>> ref = { {5,6},{6,7},{7,8},{8,5} };
//	EXPECT_EQ(ref, result);
//}
//
//TEST(Element, check_face_type_1) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 1,2 };
//
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::triangle;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n3 = { 3,3 };
//	std::vector<Euclidean_Vector<2>> cell_nodes = { n1,n2,n3 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.check_face_type(cell_element);
//	const auto ref = FaceType::inward_face;
//
//	EXPECT_EQ(result, ref);
//}
//TEST(Element, check_face_type_2) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n2,n1 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 2,1 };
//
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::triangle;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n3 = { 3,3 };
//	std::vector<Euclidean_Vector<2>> cell_nodes = { n1,n2,n3 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.check_face_type(cell_element);
//	const auto ref = FaceType::outward_face;
//
//	EXPECT_EQ(result, ref);
//}
//TEST(Element, check_face_type_3) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1,2 };
//	const Euclidean_Vector n2 = { 2,3,3 };
//	const Euclidean_Vector n3 = { 2,2,1 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 0,1,2 };
//
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::tetrahedral;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n4 = { 2,3,1 };
//	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3,n4 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 0,1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.check_face_type(cell_element);
//	const auto ref = FaceType::inward_face;
//
//	EXPECT_EQ(result, ref);
//}
//TEST(Element, check_face_type_4) {
///
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1,2 };
//	const Euclidean_Vector n2 = { 2,3,3 };
//	const Euclidean_Vector n3 = { 2,2,1 };
//	std::vector<Euclidean_Vector> nodes = { n2,n1,n3 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 1,0,2 };
//
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::tetrahedral;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n4 = { 2,3,1 };
//	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3,n4 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 0,1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.check_face_type(cell_element);
//	const auto ref = FaceType::outward_face;
//
//	EXPECT_EQ(result, ref);
//}
//
//
//TEST(Element, normalized_normal_vector_1) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 1,2 };
//	
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::triangle;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n3 = { 3,3 };
//	std::vector<Euclidean_Vector<2>> cell_nodes = { n1,n2,n3 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.normalized_normal_vector(cell_element, face_element.geometry_.center_node());
//
//	const Euclidean_Vector ref = { 0,-1 };
//
//	EXPECT_EQ(result, ref);
//}
//TEST(Element, normalized_normal_vector_2) {
///
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n2,n1 };
//
//	const auto element_type = ElementType::inner_face;
//	Geometry face_geometry(ref_geometry, std::move(nodes));
//	std::vector<uint> face_node_indexes = { 2,1 };
//
//	Element face_element(element_type, std::move(face_geometry), std::move(face_node_indexes));
//
//	const Figure cell_fig = Figure::triangle;
//	const ReferenceGeometry cell_ref_geometry(cell_fig, fig_order);
//
//	const Euclidean_Vector n3 = { 3,3 };
//	std::vector<Euclidean_Vector<2>> cell_nodes = { n1,n2,n3 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(cell_ref_geometry, std::move(cell_nodes));
//	std::vector<uint> cell_node_indexes = { 1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_geometry), std::move(cell_node_indexes));
//
//	const auto result = face_element.normalized_normal_vector(cell_element, face_element.geometry_.center_node());
//
//	const Euclidean_Vector ref = { 0,-1 };
//
//	EXPECT_EQ(result, ref);
//}
//
//TEST(ms, is_circular_permuation_1) {
//	std::vector<ushort> v1 = { 1,2,3 };
//	std::vector<ushort> v2 = { 3,1,2 };
//	
//	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
//}
//TEST(ms, is_circular_permuation_2) {
//	std::vector<ushort> v1 = { 1,2,3 };
//	std::vector<ushort> v2 = { 3,2,1 };
//
//	EXPECT_FALSE(ms::is_circular_permutation(v1, v2));
//}
//TEST(ms, is_circular_permuation_3) {
//	std::vector<ushort> v1 = { 1,2,3,4 };
//	std::vector<ushort> v2 = { 2,3,4,1 };
//
//	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
//}
//TEST(ms, is_circular_permuation_4) {
//	std::vector<ushort> v1 = { 1,2,3,4,5,6 };
//	std::vector<ushort> v2 = { 3,4,5 };
//
//	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
//}
//TEST(ms, is_circular_permuation_5) {
//	std::vector<ushort> v1 = { 1,2,3,4 };
//	std::vector<ushort> v2 = { 4,1,2,3 };
//
//	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
//}
//
//TEST(ms, contains_1) {
//	std::vector<int> v1 = { 1,2,3,4,5 };
//
//	EXPECT_TRUE(ms::contains(v1, 3));
//}
//TEST(ms, contains_2) {
//	std::vector<int> v1 = { 1,2,3,4,5 };
//
//	EXPECT_FALSE(ms::contains(v1, 6));
//}
//
//TEST(ms, has_intersection_1) {
//	std::vector<int> v1 = { 1,2,3,4,5 };
//	std::vector<int> v2 = { 66,7,8,9,1 };
//	EXPECT_TRUE(ms::has_intersection(v1, v2));
//}
//TEST(ms, has_intersection_2) {
//	std::vector<int> v1 = { 1,2,3,4,5 };
//	std::vector<int> v2 = { 66,7,8,9 };
//	EXPECT_FALSE(ms::has_intersection(v1, v2));
//}
//TEST(ms, has_intersection_3) {
//	std::vector<int> v1 = { 1,2,3,4,5 };
//	std::vector<int> v2;
//	EXPECT_FALSE(ms::has_intersection(v1, v2));
//}

//TEST(Geometry, periodic_match_2) {
//	const Figure fig = Figure::line;
//	const order_t fig_order = 1;
//	const ReferenceGeometryspace_dimension> ref_geometry(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,0 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//	
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	const Euclidean_Vector n3 = { 8,0 };
//	const Euclidean_Vector n4 = { 10,1 };
//	std::vector<Euclidean_Vector<2>> nodes2 = { n3,n4 };
//	
//
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	constexpr size_t axis_tag = 1;
//	EXPECT_FALSE(geometry.is_periodic_pair(geometry2, axis_tag));
//}
//TEST(Geometry, periodic_match_3) {
//	const Figure fig = Figure::line;
//	const order_t fig_order = 1;
//	const ReferenceGeometryspace_dimension> ref_geometry(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,0 };
//	const Euclidean_Vector n2 = { 3,1 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2 };
//	
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//
//	const Euclidean_Vector n3 = { 1,5 };
//	const Euclidean_Vector n4 = { 3,6 };
//	std::vector<Euclidean_Vector<2>> nodes2 = { n3,n4 };
//	
//
//	Geometry geometry2(ref_geometry, std::move(nodes2));
//
//	constexpr size_t axis_tag = 1;
//	EXPECT_TRUE(geometry.is_periodic_pair(geometry2, axis_tag));
//}
//
//TEST(Geometry, faces_geometry_1) {
//	const Figure fig = Figure::triangle;
//	const order_t fig_order = 1;
//	const ReferenceGeometryspace_dimension> ref_geometry(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3 };
//	
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.face_geometries();
//
//	const Figure f_fig = Figure::line;
//
//	const Euclidean_Vector f1_n1 = { 1,1 };
//	const Euclidean_Vector f1_n2 = { 2,1 };
//	std::vector<Euclidean_Vector<2>> f1_nodes = { f1_n1,f1_n2 };
//	std::vector<size_t> f1_indexes = { 1,2 };
//	Geometry f1_geometry(f_fig, fig_order, std::move(f1_nodes), std::move(f1_indexes));
//
//	const Euclidean_Vector f2_n1 = { 2,1 };
//	const Euclidean_Vector f2_n2 = { 4,2 };
//	std::vector<Euclidean_Vector<2>> f2_nodes = { f2_n1,f2_n2 };
//	std::vector<size_t> f2_indexes = { 1,2 };
//	Geometry f2_geometry(f_fig, fig_order, std::move(f2_nodes), std::move(f2_indexes));
//
//	const Euclidean_Vector f3_n1 = { 4,2 };
//	const Euclidean_Vector f3_n2 = { 1,1 };
//	std::vector<Euclidean_Vector<2>> f3_nodes = { f3_n1,f3_n2 };
//	std::vector<size_t> f3_indexes = { 1,2 };
//	Geometry f3_geometry(f_fig, fig_order, std::move(f3_nodes), std::move(f3_indexes));
//
//
//	const std::vector<Geometry<2>> ref = { f1_geometry,f2_geometry,f3_geometry };
//	EXPECT_EQ(result, ref);
//}
//
//TEST(Geometry, vertex_node_indexes_1) {
//	const Figure fig = Figure::quadrilateral;
//	const order_t fig_order = 1;
//	const ReferenceGeometryspace_dimension> ref_geometry(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector<2>> nodes = { n1,n2,n3,n4 };
//	
//
//	Geometry geometry(ref_geometry, std::move(nodes));
//	const auto result = geometry.vertex_node_indexes();
//
//	std::vector<size_t> ref = { 1,2,3,4 };
//	EXPECT_EQ(result, ref);
//}









//TEST(ReferenceGeometry, mapping_function_1) {
//
//	const Figure fig = Figure::line;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2 };
//	Euclidean_Vector p2 = { 4,2 };
//	std::vector<Euclidean_Vector> pv = { p1,p2 };
//
//	const auto result = ref_geo->get_mapping_monomial_vector_function(); .mapping_function(pv);
//
//	Euclidean_Vector rp1 = { -1,0 };
//	Euclidean_Vector rp2 = { 1,0 };
//
//	EXPECT_EQ(p1, result(rp1));
//	EXPECT_EQ(p2, result(rp2));
//}
//TEST(ReferenceGeometry, mapping_function_2) {
//
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2 };
//	Euclidean_Vector p2 = { 4,2 };
//	Euclidean_Vector p3 = { 2,3 };
//	std::vector<Euclidean_Vector> pv = { p1,p2,p3 };
//
//	const auto result = ref_geometry.mapping_function(pv);
//	//std::cout << "\n" << result << "\n";
//
//	Euclidean_Vector rp1 = { -1,-1 };
//	Euclidean_Vector rp2 = { 1,-1 };
//	Euclidean_Vector rp3 = { -1,1 };
//
//	EXPECT_EQ(p1, result(rp1));
//	EXPECT_EQ(p2, result(rp2));
//	EXPECT_EQ(p3, result(rp3));
//}
//TEST(ReferenceGeometry, mapping_function_3) {
//
//	const Figure fig = Figure::quadrilateral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2 };
//	Euclidean_Vector p2 = { 3,1 };
//	Euclidean_Vector p3 = { 4,1 };
//	Euclidean_Vector p4 = { 1,3 };
//	std::vector<Euclidean_Vector> pv = { p1,p2,p3,p4 };
//
//	const auto result = ref_geometry.mapping_function(pv);
//
//	Euclidean_Vector rp1 = { -1,-1 };
//	Euclidean_Vector rp2 = { 1,-1 };
//	Euclidean_Vector rp3 = { 1,1 };
//	Euclidean_Vector rp4 = { -1,1 };
//
//	EXPECT_EQ(p1, result(rp1));
//	EXPECT_EQ(p2, result(rp2));
//	EXPECT_EQ(p3, result(rp3));
//	EXPECT_EQ(p4, result(rp4));
//}
//TEST(ReferenceGeometry, mapping_function_4) {
///
//	const auto fig = Figure::tetrahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2,4 };
//	Euclidean_Vector p2 = { 3,1,5 };
//	Euclidean_Vector p3 = { 4,1,9 };
//	Euclidean_Vector p4 = { 1,3,2 };
//	std::vector<Euclidean_Vector> pv = { p1,p2,p3,p4 };
//
//	const auto result = ref_geometry.mapping_function(pv);
//		
//	Euclidean_Vector rp1 = { -1,-1,-1 };
//	Euclidean_Vector rp2 = { 1,-1,-1 };
//	Euclidean_Vector rp3 = { -1,1,-1 };
//	Euclidean_Vector rp4 = { -1,-1,1 };
//
//	EXPECT_EQ(p1, result(rp1));
//	EXPECT_EQ(p2, result(rp2));
//	EXPECT_EQ(p3, result(rp3));
//	EXPECT_EQ(p4, result(rp4));
//}
//TEST(ReferenceGeometry, mapping_function_5) {
///
//	const auto fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	Euclidean_Vector p1 = { 1,2,1 };
//	Euclidean_Vector p2 = { 3,1,2 };
//	Euclidean_Vector p3 = { 4,1,3 };
//	Euclidean_Vector p4 = { 1,3,1 };
//	Euclidean_Vector p5 = { 1,2,4 };
//	Euclidean_Vector p6 = { 3,1,7 };
//	Euclidean_Vector p7 = { 4,1,6 };
//	Euclidean_Vector p8 = { 1,3,8 };
//	std::vector<Euclidean_Vector> pv = { p1,p2,p3,p4,p5,p6,p7,p8 };
//
//	const auto result = ref_geometry.mapping_function(pv);
//
//	Euclidean_Vector rp1 = { -1,-1,-1 };
//	Euclidean_Vector rp2 = { 1,-1,-1 };
//	Euclidean_Vector rp3 = { 1,1,-1 };
//	Euclidean_Vector rp4 = { -1,1,-1 };
//	Euclidean_Vector rp5 = { -1,-1,1 };
//	Euclidean_Vector rp6 = { 1,-1,1 };
//	Euclidean_Vector rp7 = { 1,1,1 };
//	Euclidean_Vector rp8 = { -1,1,1 };
//
//	EXPECT_EQ(p1, result(rp1));
//	EXPECT_EQ(p2, result(rp2));
//	EXPECT_EQ(p3, result(rp3));
//	EXPECT_EQ(p4, result(rp4));
//	EXPECT_EQ(p5, result(rp5));
//	EXPECT_EQ(p6, result(rp6));
//	EXPECT_EQ(p7, result(rp7));
//	EXPECT_EQ(p8, result(rp8));
//}