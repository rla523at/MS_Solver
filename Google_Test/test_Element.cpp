#pragma once
#include "gtest/gtest.h"
#include "../MS_Solver/INC/Element.h"


TEST(Element, vertex_node_indexes_1) 
{
	const Figure fig = Figure::quadrilateral;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 2,1 };
	const Euclidean_Vector n3 = { 4,2 };
	const Euclidean_Vector n4 = { 1,2 };
	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };
	Geometry geometry(std::move(ref_geo), std::move(nodes));

	ElementType element_type = ElementType::cell;
	std::vector<uint> indexes = { 5,6,7,8 };

	Element element(element_type, std::move(indexes), std::move(geometry));
	const auto result = element.vertex_node_indexes();

	const std::vector<uint> ref = { 5,6,7,8 };
	EXPECT_EQ(ref, result);
}

//TEST(Element, set_of_face_node_indexes) 
//{
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
//	ElementType element_type = ElementType::cell;
//	std::vector<uint> indexes = { 5,6,7,8 };
//
//	Element element(element_type, std::move(indexes), std::move(geometry));
//	const auto result = element.set_of_face_node_indexes();
//
//	const std::vector<std::vector<uint>> ref = { {5,6},{6,7},{7,8},{8,5} };
//	EXPECT_EQ(ref, result);
//}
//TEST(Element, set_of_face_vertex_node_indexes) 
//{
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
//	ElementType element_type = ElementType::cell;
//	std::vector<uint> indexes = { 5,6,7,8 };
//
//	Element element(element_type, std::move(indexes), std::move(geometry));
//
//	const auto result = element.set_of_face_vertex_node_indexes();
//
//	const std::vector<std::vector<uint>> ref = { {5,6},{6,7},{7,8},{8,5} };
//	EXPECT_EQ(ref, result);
//}

TEST(Element, outward_normalized_normal_vector_1) 
{
	//make face element
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };

	const auto element_type = ElementType::face;
	Geometry face_geometry(std::move(ref_geo), std::move(nodes));
	std::vector<uint> face_node_indexes = { 1,2 };

	Element face_element(element_type, std::move(face_node_indexes), std::move(face_geometry));

	//make cell element
	const Figure cell_fig = Figure::triangle;
	auto cell_ref_geo = Reference_Geometry_Factory::make(cell_fig, fig_order);

	const Euclidean_Vector n3 = { 3,3 };
	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3 };

	const auto cell_element_type = ElementType::cell;
	Geometry cell_geometry(std::move(cell_ref_geo), std::move(cell_nodes));
	std::vector<uint> cell_node_indexes = { 1,2,3 };

	Element cell_element(cell_element_type, std::move(cell_node_indexes), std::move(cell_geometry));

	const auto result = face_element.outward_normalized_normal_vector(cell_element, face_element.center_point());
	Euclidean_Vector ref = { 0,-1 };
	EXPECT_EQ(result, ref);
}
TEST(Element, outward_normalized_normal_vector_2) 
{
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,1 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n2,n1 };

	const auto element_type = ElementType::face;
	Geometry face_geometry(std::move(ref_geo), std::move(nodes));
	std::vector<uint> face_node_indexes = { 2,1 };

	Element face_element(element_type, std::move(face_node_indexes), std::move(face_geometry));


	const Figure cell_fig = Figure::triangle;
	auto cell_ref_geo = Reference_Geometry_Factory::make(cell_fig, fig_order);


	const Euclidean_Vector n3 = { 3,3 };
	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3 };

	const auto cell_element_type = ElementType::cell;
	Geometry cell_geometry(std::move(cell_ref_geo), std::move(cell_nodes));

	std::vector<uint> cell_node_indexes = { 1,2,3 };

	Element cell_element(cell_element_type, std::move(cell_node_indexes), std::move(cell_geometry));

	const auto result = face_element.outward_normalized_normal_vector(cell_element, face_element.center_point());
	Euclidean_Vector ref = { 0,-1 };
	EXPECT_EQ(result, ref);
}
//TEST(Element, check_face_type_3) 
//{
//	//make face element
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1,2 };
//	const Euclidean_Vector n2 = { 2,3,3 };
//	const Euclidean_Vector n3 = { 2,2,1 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3 };
//
//	Geometry face_geometry(std::move(ref_geo), std::move(nodes));
//
//	const auto element_type = ElementType::face;
//	std::vector<uint> face_node_indexes = { 0,1,2 };
//	Element face_element(element_type, std::move(face_node_indexes), std::move(face_geometry));
//
//	//make cell element
//	const Figure cell_fig = Figure::tetrahedral;
//	auto cell_ref_geo = Reference_Geometry_Factory::make(cell_fig, fig_order);
//
//	const Euclidean_Vector n4 = { 2,3,1 };
//	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3,n4 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(std::move(cell_ref_geo), std::move(cell_nodes));
//
//	std::vector<uint> cell_node_indexes = { 0,1,2,3 };
//	Element cell_element(cell_element_type, std::move(cell_node_indexes), std::move(cell_geometry));
//
//	const auto result = face_element.outward_normalized_normal_vector(cell_element, face_element.center_node());
//	Euclidean_Vector ref = { 0.801783725737273, -0.534522483824849, -0.267261241912424 };
//	EXPECT_EQ(result, ref);
//}
//TEST(Element, check_face_type_4) {
//
//	const Figure fig = Figure::triangle;
//	const ushort fig_order = 1;
//	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1,2 };
//	const Euclidean_Vector n2 = { 2,3,3 };
//	const Euclidean_Vector n3 = { 2,2,1 };
//	std::vector<Euclidean_Vector> nodes = { n2,n1,n3 };
//
//	const auto element_type = ElementType::face;
//
//	Geometry face_geometry(std::move(ref_geo), std::move(nodes));
//
//	std::vector<uint> face_node_indexes = { 1,0,2 };
//
//	Element face_element(element_type, std::move(face_node_indexes), std::move(face_geometry));
//
//
//	const Figure cell_fig = Figure::tetrahedral;
//	auto cell_ref_geo = Reference_Geometry_Factory::make(cell_fig, fig_order);
//
//
//	const Euclidean_Vector n4 = { 2,3,1 };
//	std::vector<Euclidean_Vector> cell_nodes = { n1,n2,n3,n4 };
//
//	const auto cell_element_type = ElementType::cell;
//	Geometry cell_geometry(std::move(cell_ref_geo), std::move(cell_nodes));
//
//	std::vector<uint> cell_node_indexes = { 0,1,2,3 };
//
//	Element cell_element(cell_element_type, std::move(cell_node_indexes), std::move(cell_geometry));
//
//
//	const auto result = face_element.check_face_type(cell_element);
//	const auto ref = FaceType::outward_face;
//
//	EXPECT_EQ(result, ref);
//}

TEST(Element, find_periodic_matched_node_indexes_1)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,0 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 1,2 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 8,0 };
	const Euclidean_Vector n4 = { 10,1 };
	std::vector<Euclidean_Vector> nodes2 = { n3,n4 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 3,4 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { 3,4 };

	EXPECT_EQ(result, ref);
}
TEST(Element, find_periodic_matched_node_indexes_2)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,0 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 1,2 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 8,0 };
	const Euclidean_Vector n4 = { 10,1 };
	std::vector<Euclidean_Vector> nodes2 = { n4,n3 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 4,3 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { 3,4 };

	EXPECT_EQ(result, ref);
}
TEST(Element, find_periodic_matched_node_indexes_3)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,0 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 1,2 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 1,5 };
	const Euclidean_Vector n4 = { 3,6 };
	std::vector<Euclidean_Vector> nodes2 = { n4,n3 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 4,3 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { 3,4 };

	EXPECT_EQ(result, ref);
}
TEST(Element, find_periodic_matched_node_indexes_4)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1, 0 };
	const Euclidean_Vector n2 = { 1, 1.0 / 3.0 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 3,7 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 0, 1.0 / 3.0 };
	const Euclidean_Vector n4 = { 0, 0 };
	std::vector<Euclidean_Vector> nodes2 = { n3,n4 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 4,0 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { 0,4 };

	EXPECT_EQ(result, ref);
}
TEST(Element, find_periodic_matched_node_indexes_5)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1, 0, 0 };
	const Euclidean_Vector n2 = { 1, 1.0 / 3.0, 0 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 3,7 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 0, 1.0 / 3.0, 0 };
	const Euclidean_Vector n4 = { 0, 0, 0 };
	std::vector<Euclidean_Vector> nodes2 = { n3,n4 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 4,0 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { }; // 두 element는 z=0 plane에 같이 있음으로 periodic pair가 될 수 없다.

	EXPECT_EQ(result, ref);
}
TEST(Element, find_periodic_matched_node_indexes_6)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1, 0, 0 };
	const Euclidean_Vector n2 = { 1, 1.0 / 3.0, 0 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 3,7 };
	Element element1(element_type, std::move(node_indexes), std::move(geometry1));

	//make element 2
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n3 = { 0, 1.0 / 3.0, 1 };
	const Euclidean_Vector n4 = { 0, 0, 1 };
	std::vector<Euclidean_Vector> nodes2 = { n3,n4 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 4,0 };
	Element element2(element_type, std::move(node_indexes2), std::move(geometry2));

	const auto result = element2.find_periodic_matched_node_indexes(element1);
	std::vector<uint> ref = { }; // 두 element는 axis translation이 되지 않아 periodic pair가 될 수 없다.

	EXPECT_EQ(result, ref);
}

TEST(Element, rearrange_node_indexes_1)
{
	//make element 1
	const Figure fig = Figure::line;
	const ushort fig_order = 1;
	auto ref_geo = Reference_Geometry_Factory::make(fig, fig_order);

	const Euclidean_Vector n1 = { 1,0 };
	const Euclidean_Vector n2 = { 3,1 };
	std::vector<Euclidean_Vector> nodes = { n1,n2 };
	Geometry geometry1(std::move(ref_geo), std::move(nodes));

	const auto element_type = ElementType::periodic;
	std::vector<uint> node_indexes = { 1,2 };
	Element element(element_type, std::move(node_indexes), std::move(geometry1));

	element.rearrange_node_indexes({ 2,1 });

	//make ref element
	auto ref_geo2 = Reference_Geometry_Factory::make(fig, fig_order);
	std::vector<Euclidean_Vector> nodes2 = { n2,n1 };
	Geometry geometry2(std::move(ref_geo2), std::move(nodes2));

	std::vector<uint> node_indexes2 = { 2,1 };
	Element ref_element(element_type, std::move(node_indexes2), std::move(geometry2));

	EXPECT_EQ(element, ref_element);
}

TEST(ms, is_circular_permuation_1) 
{
	std::vector<ushort> v1 = { 1,2,3 };
	std::vector<ushort> v2 = { 3,1,2 };
	
	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
}
TEST(ms, is_circular_permuation_2) 
{
	std::vector<ushort> v1 = { 1,2,3 };
	std::vector<ushort> v2 = { 3,2,1 };

	EXPECT_FALSE(ms::is_circular_permutation(v1, v2));
}
TEST(ms, is_circular_permuation_3) 
{
	std::vector<ushort> v1 = { 1,2,3,4 };
	std::vector<ushort> v2 = { 2,3,4,1 };

	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
}
TEST(ms, is_circular_permuation_4) 
{
	std::vector<ushort> v1 = { 1,2,3,4,5,6 };
	std::vector<ushort> v2 = { 3,4,5 };

	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
}
TEST(ms, is_circular_permuation_5) 
{
	std::vector<ushort> v1 = { 1,2,3,4 };
	std::vector<ushort> v2 = { 4,1,2,3 };

	EXPECT_TRUE(ms::is_circular_permutation(v1, v2));
}

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





//TEST(Geometry, vertex_node_indexes_1) {
//	const Figure fig = Figure::quadrilateral;
//	const order_t fig_order = 1;
//	const ReferenceGeometryspace_dimension> ref_geometry(fig, fig_order);
//
//	const Euclidean_Vector n1 = { 1,1 };
//	const Euclidean_Vector n2 = { 2,1 };
//	const Euclidean_Vector n3 = { 4,2 };
//	const Euclidean_Vector n4 = { 1,2 };
//	std::vector<Euclidean_Vector> nodes = { n1,n2,n3,n4 };
//	
//
//	Geometry geometry(std::move(ref_geo), std::move(nodes));
//	const auto result = geometry.vertex_node_indexes();
//
//	std::vector<size_t> ref = { 1,2,3,4 };
//	EXPECT_EQ(result, ref);
//}
//
//






