#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Reference_Geometry_Impl.h"

TEST(Reference_Geometry, mapping_monomial_vector_function_1) 
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_monomial_vector_function();

	const Polynomial r("x0");
	const Vector_Function<Polynomial> ref = { 1, r };
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, mapping_monomial_vector_function_2) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_monomial_vector_function();

	const Polynomial r("x0");
	const Polynomial s("x1");

	const Vector_Function<Polynomial> ref = { 1, r, s };
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, mapping_monomial_vector_function_3) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_monomial_vector_function();

	const Polynomial r("x0");
	const Polynomial s("x1");

	const Vector_Function<Polynomial> ref = { 1, r, r * s, s };
	EXPECT_EQ(ref, result);
}
//TEST(Reference_Geometry, mapping_monomial_vector_function_4) 
//{
//	const Figure fig = Figure::tetrahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//	const auto& result = ref_geo->get_mapping_monomial_vector_function();
//
//	const Polynomial r("x0");
//	const Polynomial s("x1");
//	const Polynomial t("x2");
//
//	const Vector_Function<Polynomial> ref = { 1, r, s, t };
//	EXPECT_EQ(ref, result);
//}
//TEST(Reference_Geometry, mapping_monomial_vector_function_5) 
//{
//	const Figure fig = Figure::tetrahedral;
//	const ushort fig_order = 2;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//	const auto& result = ref_geo->get_mapping_monomial_vector_function();
//
//	const Polynomial r("x0");
//	const Polynomial s("x1");
//	const Polynomial t("x2");
//
//	const Vector_Function<Polynomial> ref = { 1, r, s, t, r ^ 2, r * s , r * t, s ^ 2, s * t, t ^ 2 };
//	EXPECT_EQ(ref, result);
//}
//TEST(Reference_Geometry, mapping_monomial_vector_function_6) 
//{
//	const Figure fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//	const auto& result = ref_geo->get_mapping_monomial_vector_function();
//
//	const Polynomial r("x0");
//	const Polynomial s("x1");
//	const Polynomial t("x2");
//
//	const Vector_Function<Polynomial> ref = { 1, r, r * s, s, t, r * t, r * s * t, s * t };
//	EXPECT_EQ(ref, result);
//}

TEST(Reference_Geometry, is_simplex1) 
{
	const auto fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	EXPECT_TRUE(ref_geo->is_simplex());
}
TEST(Reference_Geometry, is_simplex2) 
{
	const auto fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	EXPECT_FALSE(ref_geo->is_simplex());
}
//TEST(Reference_Geometry, is_simplex3) 
//{
//	const auto fig = Figure::tetrahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	EXPECT_TRUE(ref_geo->is_simplex());
//}
//TEST(Reference_Geometry, is_simplex4) 
//{
//	const auto fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	EXPECT_FALSE(ref_geo->is_simplex());
//}

TEST(Reference_Geometry, get_mapping_nodes1) 
{
	const auto fig = Figure::line;
	const auto fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_nodes();

	const std::vector<Euclidean_Vector> ref = { { -1,0 }, { 1,0 } };
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, get_mapping_nodes2) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_nodes();

	const std::vector<Euclidean_Vector> ref = { { -1, -1 }, { 1, -1 }, { -1, 1 } };
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, get_mapping_nodes3) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_mapping_nodes();

	const std::vector<Euclidean_Vector> ref = { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, inverse_mapping_monomial_matrix_1) 
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_inverse_mapping_monomial_matrix();

	const Matrix ref(2, 2, { 0.5,-0.5,0.5,0.5 });
	EXPECT_EQ(ref, result);
}
TEST(Reference_Geometry, inverse_mapping_monomial_matrix_2) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_inverse_mapping_monomial_matrix();

	const Matrix ref(3, 3, { 0, -0.5, -0.5,0.5,  0.5,    0,0.5,    0,  0.5 });
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, inverse_mapping_monomial_matrix_3) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);
	const auto& result = ref_geo->get_inverse_mapping_monomial_matrix();

	const Matrix ref(4, 4, { 0.25, -0.25,  0.25, -0.25,  0.25,  0.25, -0.25, -0.25,  0.25,  0.25,  0.25,  0.25,  0.25, -0.25, -0.25,  0.25 });
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_quadrature_rule_1) 
{
	const auto fig = Figure::line;
	const auto fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr size_t integrand_order = 0;
	const auto result = ref_geo->get_quadrature_rule(integrand_order);

	Quadrature_Rule ref = { { { 0.000000000000000 } }, { 2.000000000000000 } };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_quadrature_rule_2)
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	for (ushort i = 0; i < 22; ++i)
	{
		const auto ref_quad_rule = ref_geo->get_quadrature_rule(i);

		double sum = 0.0;
		for (const auto& weight : ref_quad_rule.weights)
			sum += weight;

		const auto ref = 2.0;
		//EXPECT_DOUBLE_EQ(sum, ref);	//round off error
		EXPECT_NEAR(sum, ref, 9.0E-15);
	}
}
TEST(Reference_Geometry, get_quadrature_rule_3) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr size_t integrand_order = 11;
	const auto result = ref_geo->get_quadrature_rule(integrand_order);

	Quadrature_Rule ref = { { { -0.9504029146358142, -0.949107912342758 }, { -0.9556849211223275, -0.741531185599394 }, { -0.9642268026617964, -0.405845151377397 }, { -0.974553956171379, 0 }, { -0.9848811096809615, 0.405845151377397 }, { -0.9934229912204304, 0.741531185599394 }, { -0.9987049977069438, 0.949107912342758 }, { -0.7481081943789636, -0.949107912342758 }, { -0.7749342496082214, -0.741531185599394 }, { -0.8183164352463219, -0.405845151377397 }, { -0.870765592799697, 0 }, { -0.9232147503530721, 0.405845151377397 }, { -0.9665969359911726, 0.741531185599394 }, { -0.9934229912204304, 0.949107912342758 }, { -0.4209640416964355, -0.949107912342758 }, { -0.4826304010243249, -0.741531185599394 }, { -0.5823551434482712, -0.405845151377397 }, { -0.7029225756886985, 0 }, { -0.8234900079291259, 0.405845151377397 }, { -0.9232147503530722, 0.741531185599394 }, { -0.9848811096809615, 0.949107912342758 }, { -0.02544604382862098, -0.949107912342758 }, { -0.129234407200303, -0.741531185599394 }, { -0.2970774243113015, -0.405845151377397 }, { -0.5, 0 }, { -0.7029225756886985, 0.405845151377397 }, { -0.870765592799697, 0.741531185599394 }, { -0.974553956171379, 0.949107912342758 }, { 0.3700719540391935, -0.949107912342758 }, { 0.2241615866237189, -0.741531185599394 }, { -0.01179970517433182, -0.405845151377397 }, { -0.2970774243113015, 0 }, { -0.5823551434482711, 0.405845151377397 }, { -0.8183164352463218, 0.741531185599394 }, { -0.9642268026617964, 0.949107912342758 }, { 0.6972161067217216, -0.949107912342758 }, { 0.5164654352076156, -0.741531185599394 }, { 0.2241615866237189, -0.405845151377397 }, { -0.129234407200303, 0 }, { -0.4826304010243249, 0.405845151377397 }, { -0.7749342496082215, 0.741531185599394 }, { -0.9556849211223276, 0.949107912342758 }, { 0.8995108269785723, -0.949107912342758 }, { 0.6972161067217215, -0.741531185599394 }, { 0.3700719540391935, -0.405845151377397 }, { -0.02544604382862098, 0 }, { -0.4209640416964355, 0.405845151377397 }, { -0.7481081943789636, 0.741531185599394 }, { -0.9504029146358142, 0.949107912342758 } }, { 0.01633971902233046, 0.03153707751100931, 0.03475337161903316, 0.02705971537896383, 0.01468787955288011, 0.004680565643230263, 0.0004266374414229519, 0.03529604741916743, 0.06812443847836634, 0.07507207749196593, 0.05845271854796313, 0.03172784626693878, 0.01011066754980336, 0.0009215957350721348, 0.04818316692765089, 0.09299769892288511, 0.1024820257759689, 0.07979468810555948, 0.04331216169277281, 0.0138022248360196, 0.001258084244262363, 0.05274230535088141, 0.1017972322343419, 0.1121789753788724, 0.08734493960849631, 0.04741040083224651, 0.01510820486158434, 0.001377125407046246, 0.04818316692765089, 0.09299769892288511, 0.1024820257759689, 0.07979468810555948, 0.04331216169277281, 0.0138022248360196, 0.001258084244262363, 0.03529604741916743, 0.06812443847836634, 0.07507207749196593, 0.05845271854796313, 0.03172784626693878, 0.01011066754980336, 0.0009215957350721348, 0.01633971902233046, 0.03153707751100931, 0.03475337161903316, 0.02705971537896383, 0.01468787955288011, 0.004680565643230263, 0.0004266374414229519 } };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_quadrature_rule_4) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	for (ushort i = 0; i < 21; ++i) 
{
		const auto ref_quad_rule = ref_geo->get_quadrature_rule(i);

		double sum = 0.0;
		for (const auto& weight : ref_quad_rule.weights)
			sum += weight;

		const auto ref = 2.0;
		//EXPECT_DOUBLE_EQ(sum, ref);	//round off error
		EXPECT_NEAR(sum, ref, 9.0E-15);
	}
}
TEST(Reference_Geometry, get_quadrature_rule_5) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr size_t integrand_order = 11;
	const auto result = ref_geo->get_quadrature_rule(integrand_order);

	Quadrature_Rule ref = { { { -0.9324695142031521, -0.9324695142031521 }, { -0.9324695142031521, -0.661209386466264 }, { -0.9324695142031521, -0.238619186083197 }, { -0.9324695142031521, 0.238619186083197 }, { -0.9324695142031521, 0.661209386466264 }, { -0.9324695142031521, 0.9324695142031521 }, { -0.661209386466264, -0.9324695142031521 }, { -0.661209386466264, -0.661209386466264 }, { -0.661209386466264, -0.238619186083197 }, { -0.661209386466264, 0.238619186083197 }, { -0.661209386466264, 0.661209386466264 }, { -0.661209386466264, 0.9324695142031521 }, { -0.238619186083197, -0.9324695142031521 }, { -0.238619186083197, -0.661209386466264 }, { -0.238619186083197, -0.238619186083197 }, { -0.238619186083197, 0.238619186083197 }, { -0.238619186083197, 0.661209386466264 }, { -0.238619186083197, 0.9324695142031521 }, { 0.238619186083197, -0.9324695142031521 }, { 0.238619186083197, -0.661209386466264 }, { 0.238619186083197, -0.238619186083197 }, { 0.238619186083197, 0.238619186083197 }, { 0.238619186083197, 0.661209386466264 }, { 0.238619186083197, 0.9324695142031521 }, { 0.661209386466264, -0.9324695142031521 }, { 0.661209386466264, -0.661209386466264 }, { 0.661209386466264, -0.238619186083197 }, { 0.661209386466264, 0.238619186083197 }, { 0.661209386466264, 0.661209386466264 }, { 0.661209386466264, 0.9324695142031521 }, { 0.9324695142031521, -0.9324695142031521 }, { 0.9324695142031521, -0.661209386466264 }, { 0.9324695142031521, -0.238619186083197 }, { 0.9324695142031521, 0.238619186083197 }, { 0.9324695142031521, 0.661209386466264 }, { 0.9324695142031521, 0.9324695142031521 } }, { 0.02935208168898062, 0.06180729337238363, 0.08016511731780691, 0.08016511731780691, 0.06180729337238363, 0.02935208168898062, 0.06180729337238363, 0.1301489125881677, 0.168805367087588, 0.168805367087588, 0.1301489125881677, 0.06180729337238363, 0.08016511731780691, 0.168805367087588, 0.2189434501672965, 0.2189434501672965, 0.168805367087588, 0.08016511731780691, 0.08016511731780691, 0.168805367087588, 0.2189434501672965, 0.2189434501672965, 0.168805367087588, 0.08016511731780691, 0.06180729337238363, 0.1301489125881677, 0.168805367087588, 0.168805367087588, 0.1301489125881677, 0.06180729337238363, 0.02935208168898062, 0.06180729337238363, 0.08016511731780691, 0.08016511731780691, 0.06180729337238363, 0.02935208168898062 } };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_quadrature_rule_6) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	for (ushort i = 0; i < 21; ++i) 
{
		const auto ref_quad_rule = ref_geo->get_quadrature_rule(i);

		double sum = 0.0;
		for (const auto& weight : ref_quad_rule.weights)
			sum += weight;

		const auto ref = 4.0;
		//EXPECT_DOUBLE_EQ(sum, ref);	//round off error
		EXPECT_NEAR(sum, ref, 1.0E-13);
	}
}
//TEST(Reference_Geometry, quadrature_rule7) 
//{
//
//	const auto fig = Figure::tetrahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	for (ushort i = 0; i < 11; ++i) 
//{
//		const auto ref_quad_rule = ref_geometry.quadrature_rule(i);
//
//		double sum = 0.0;
//		for (const auto& weight : ref_quad_rule.weights)
//			sum += weight;
//
//		const auto ref = 8.0/6.0;
//		//EXPECT_DOUBLE_EQ(sum, ref);	//round off error
//		EXPECT_NEAR(sum, ref, 1.0E-13);
//	}
//}
//TEST(Reference_Geometry, quadrature_rule8) 
//{
//
//	const auto fig = Figure::hexahedral;
//	const ushort fig_order = 1;
//	const auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	for (ushort i = 0; i < 21; ++i) 
//{
//		const auto ref_quad_rule = ref_geometry.quadrature_rule(i);
//
//		double sum = 0.0;
//		for (const auto& weight : ref_quad_rule.weights)
//			sum += weight;
//
//		const auto ref = 8.0;
//		//EXPECT_DOUBLE_EQ(sum, ref);	//round off error
//		EXPECT_NEAR(sum, ref, 1.0E-13);
//	}
//}

TEST(Reference_Geometry, get_connectivities1) 
{
	const auto fig = Figure::triangle;
	const auto fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr ushort post_order = 0;
	const auto& result = ref_geo->get_connectivities(post_order);

	const std::vector<std::vector<uint>> ref = { {0,1,2} };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_connectivities2) 
{
	const Figure fig = Figure::triangle;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr ushort post_order = 1;
	const auto& result = ref_geo->get_connectivities(post_order);

	const std::vector<std::vector<uint>> ref = { {0,1,3},{1,4,3},{1,2,4},{3,4,5} };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_connectivities3) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr ushort post_order = 0;
	const auto& result = ref_geo->get_connectivities(post_order);

	const std::vector<std::vector<uint>> ref = { {0,1,2},{1,3,2} };
	EXPECT_EQ(result, ref);
}
TEST(Reference_Geometry, get_connectivities4) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	const auto ref_geo = Reference_Geometry_Factory::make_unique(fig, fig_order);

	constexpr ushort post_order = 1;
	const auto& result = ref_geo->get_connectivities(post_order);

	const std::vector<std::vector<uint>> ref = { {0,1,3},{1,4,3},{1,2,4},{2,5,4},{3,4,6},{4,7,6},{4,5,7},{5,8,7} };
	EXPECT_EQ(result, ref);
}
//TEST(Reference_Geometry, scale_function_1) 
//{
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
//	const auto mapping_function = ref_geometry->mapping_function(pv);
//	const auto result = ref_geometry.scale_function(mapping_function);
//
//	const auto x = Polynomial("x0");
//	const auto y = Polynomial("x1");
//	Irrational_Function ref = (0.25 * y + 1.25) * (-0.25 * x + 0.25) - (0.25 * x + 0.25) * (-0.25 * y - 0.75);
//	EXPECT_EQ(result, ref);
//}

TEST(Reference_Geometry, face_reference_geometry_1)
{
	const auto figure = Figure::triangle;
	const auto figure_order = 1;
	const auto reference_geometry = Reference_Geometry_Factory::make_unique(figure, figure_order);

	const auto result = reference_geometry->face_reference_geometries();

	const auto line = Reference_Geometry_Factory::make_unique(Figure::line, figure_order);
	for (ushort i = 0; i < 3; ++i) 
{
		EXPECT_EQ(typeid(*result[i]), typeid(*line));
	}	
}

//TEST(Reference_Geometry, mapping_function_1) 
//{
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
//TEST(Reference_Geometry, mapping_function_2) 
//{
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
//TEST(Reference_Geometry, mapping_function_3) 
//{
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
//TEST(Reference_Geometry, mapping_function_4) 
//{
//
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
//TEST(Reference_Geometry, mapping_function_5) 
//{
//
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