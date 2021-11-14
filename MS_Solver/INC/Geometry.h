#pragma once
#include "Reference_Geometry.h"

class Geometry
{
public:
	Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes); // 의존성이 들어나도록 생성자 지우기 --> 지울필요까지 있을까 ?
	Geometry(std::unique_ptr<Reference_Geometry>&& reference_goemetry, std::vector<Euclidean_Vector>&& consisting_nodes);

public:
	Euclidean_Vector center_node(void) const;
	std::vector<Geometry> face_geometries(void) const;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
	bool is_line(void) const;
	ushort num_post_nodes(const ushort post_order) const;
	ushort num_post_elements(const ushort post_order) const;
	Euclidean_Vector normalized_normal_vector(const Euclidean_Vector& node) const;
	Vector_Function<Polynomial> orthonormal_basis_vector_function(const ushort solution_order) const;
	std::vector<Euclidean_Vector> post_nodes(const ushort post_order) const;
	std::vector<std::vector<int>> post_connectivities(const ushort post_order, const size_t connectivity_start_index) const;
	std::vector<double> projected_volume(void) const;
	std::vector<std::vector<Euclidean_Vector>> set_of_face_nodes(void) const;
	double volume(void) const;

protected:
	bool can_be_periodic_pair(const Geometry& other) const;
	ushort check_space_dimension(void) const;
	Vector_Function<Polynomial> initial_basis_vector_function(const ushort solution_order) const;
	bool is_axis_parallel_node(const Euclidean_Vector& node) const;
	bool is_on_axis(const ushort axis_tag) const;
	bool is_on_same_axis(const Geometry& other) const;
	Vector_Function<Polynomial> make_mapping_function(void) const;
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const;

protected:
	ushort space_dimension_ = 0;
	std::unique_ptr<Reference_Geometry> reference_geometry_;
	std::vector<Euclidean_Vector> nodes_;
	Vector_Function<Polynomial> mapping_function_;
	Irrational_Function scale_function_;
	mutable std::map<size_t, Quadrature_Rule> integrand_order_to_quadrature_rule_;


//public:
//	bool operator==(const Geometry & other) const;
//
//public:
//	std::vector<Space_Vector_> vertex_nodes(void) const;
//	std::vector<Geometry> sub_simplex_geometries(void) const;
//	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
//
//	template <ushort polynomial_order>
//	auto initial_basis_function(void) const;
//	template <ushort polynomial_order>
//	auto orthonormal_basis_vector_function(void) const;
//
//	//auto initial_basis_function(const ushort polynomial_order) const;
//	//auto orthonormal_basis_vector_function(const ushort polynomial_order) const;
//
//public:
//	const std::vector<Space_Vector_>& get_nodes(void) const { return nodes_; };

		//std::vector<std::vector<Space_Vector_>> set_of_face_nodes(void) const;
};


template <typename Function>
Vector_Function<Function> operator*(const Matrix& matrix, const Vector_Function<Function>& vector_function) {
	const auto [num_row, num_column] = matrix.size();

	std::vector<Function> functions(num_row);

	for (size_t i = 0; i < num_row; ++i)
		for (size_t j = 0; j < num_column; ++j)
			functions[i] += matrix.at(i, j) * vector_function[j];

	return functions;
}

namespace ms
{	
	template <typename T, typename Container>	std::vector<T> extract_by_index(const std::vector<T>& set, const Container& indexes) {
		const auto num_extracted_value = indexes.size();

		std::vector<T> extracted_values;
		extracted_values.reserve(num_extracted_value);

		for (const auto& index : indexes)
			extracted_values.push_back(set[index]);

		return extracted_values;
	}
	double integrate(const Polynomial& integrand, const Quadrature_Rule& quadrature_rule);	
	double integrate(const Polynomial& integrand, const Geometry& geometry);	
	double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry);	
	double L2_Norm(const Polynomial& function, const Geometry& geometry);
	Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& functions, const Geometry& geometry);	
}




























//Euclidean_Vector Geometry::center_node(void) const {
//	return this->mapping_function_(this->reference_geometry_.center_node());
//}
//

//std::vector<Euclidean_Vector> Geometry::post_nodes(const ushort post_order) const {
//	return this->reference_geometry_.post_nodes(this->mapping_function_, post_order);
//}
//

//std::vector<Euclidean_Vector> Geometry::vertex_nodes(void) const {
//	const auto num_vertex = this->reference_geometry_.num_vertex();
//
//	std::vector<Euclidean_Vector> vertex_nodes(this->nodes_.begin(), this->nodes_.begin() + num_vertex);
//	return vertex_nodes;
//}
//

//double Geometry::volume(void) const {
//	const auto& quadrature_rule = this->get_quadrature_rule(0);
//
//	auto volume = 0.0;
//	for (const auto weight : quadrature_rule.weights)
//		volume += weight;
//
//	return volume;
//}
//

//std::array<double, space_dimension> Geometry::projected_volume(void) const {
//	//This only work for linear mesh and convex geometry.
//
//	if constexpr (space_dimension == 2) {
//		double x_projected_volume = 0.0;
//		double y_projected_volume = 0.0;
//
//		const auto set_of_face_nodes = this->set_of_face_nodes();
//		for (const auto& face_nodes : set_of_face_nodes) {
//			const auto& start_node = face_nodes[0];
//			const auto& end_node = face_nodes[1];
//			const auto node_to_node = end_node - start_node;
//
//			x_projected_volume += std::abs(node_to_node.at(0));
//			y_projected_volume += std::abs(node_to_node.at(1));
//		}
//
//		return { 0.5 * y_projected_volume, 0.5 * x_projected_volume };
//	}
//	else if constexpr (space_dimension == 3) {
//		double yz_projected_volume = 0.0;
//		double xz_projected_volume = 0.0;
//		double xy_projected_volume = 0.0;
//
//		const auto face_geometries = this->face_geometries();
//		for (const auto& geometry : face_geometries) {
//			const auto normal_vector = geometry.normalized_normal_vector(geometry.center_node());
//
//			Euclidean_Vector yz_plane_normalized_normal_vector = { 1,0,0 };
//			Euclidean_Vector xz_plane_normalized_normal_vector = { 0,1,0 };
//			Euclidean_Vector xy_plane_normalized_normal_vector = { 0,0,1 };
//
//			const auto volume = geometry.volume();
//
//			yz_projected_volume += volume * std::abs(normal_vector.inner_product(yz_plane_normalized_normal_vector));
//			xz_projected_volume += volume * std::abs(normal_vector.inner_product(xz_plane_normalized_normal_vector));
//			xy_projected_volume += volume * std::abs(normal_vector.inner_product(xy_plane_normalized_normal_vector));
//		}
//
//		return { 0.5 * yz_projected_volume, 0.5 * xz_projected_volume, 0.5 * xy_projected_volume };
//	}
//	else {
//		throw std::runtime_error("not supproted dimension");
//		return {};
//	}
//
//
//
//}
//

//Euclidean_Vector Geometry::normalized_normal_vector(const Space_Vector_& node) const {
//	const auto normal_vector_function = this->reference_geometry_.normal_vector_function(this->mapping_function_);
//	auto normal_vector = normal_vector_function(node);
//
//	return normal_vector.be_normalize();
//}
//

//std::vector<Geometry> Geometry::face_geometries(void) const {
	//const auto face_reference_geometries = this->reference_geometry_.face_reference_geometries();
	//const auto set_of_face_node_index_orders = this->reference_geometry_.set_of_face_node_index_orders();
	//const auto num_face = set_of_face_node_index_orders.size();

	//std::vector<Geometry> face_geometries;
	//face_geometries.reserve(num_face);

	//for (size_t i = 0; i < num_face; ++i) {
	//	const auto& reference_geometry = face_reference_geometries[i];
	//	auto face_nodes = ms::extract_by_index(this->nodes_, set_of_face_node_index_orders[i]);

	//	face_geometries.push_back({ reference_geometry, std::move(face_nodes) });
	//}

	//return face_geometries;
//}
//

//std::vector<Geometry> Geometry::sub_simplex_geometries(void) const {
//	const auto sub_simplex_reference_geometries = this->reference_geometry_.sub_simplex_reference_geometries();
//	const auto set_of_sub_simplex_vnode_index_orders = this->reference_geometry_.set_of_sub_simplex_vertex_node_index_orders();
//	const auto num_sub_simplex = set_of_sub_simplex_vnode_index_orders.size();
//
//	std::vector<Geometry> sub_simplex_geometries;
//	sub_simplex_geometries.reserve(num_sub_simplex);
//
//	for (ushort i = 0; i < num_sub_simplex; ++i) {
//		const auto& reference_geometry = sub_simplex_reference_geometries[i];
//		auto sub_simplex_vnodes = ms::extract_by_index(this->nodes_, set_of_sub_simplex_vnode_index_orders[i]);
//
//		sub_simplex_geometries.push_back({ reference_geometry,std::move(sub_simplex_vnodes) });
//	}
//
//	return sub_simplex_geometries;
//}
//
//
//
//

//std::vector<std::vector<Euclidean_Vector>> Geometry::set_of_face_nodes(void) const {
//	const auto set_of_face_node_index_orders = this->reference_geometry_.set_of_face_node_index_orders();
//	const auto num_face = set_of_face_node_index_orders.size();
//
//	std::vector<std::vector<Space_Vector_>> set_of_face_nodes;
//	set_of_face_nodes.reserve(num_face);
//
//	for (size_t i = 0; i < num_face; ++i) 
//		set_of_face_nodes.push_back(ms::extract_by_index(this->nodes_, set_of_face_node_index_orders[i]));
//
//	return set_of_face_nodes;
//}
//

//

//template <ushort polynomial_order>
//auto Geometry::initial_basis_function(void) const {
//	constexpr auto num_basis = ms::combination_with_repetition(1 + space_dimension, polynomial_order);
//
//	std::array<Polynomial, num_basis> initial_basis_set = { 0 };
//
//	ushort index = 0;
//	if (space_dimension == 2) {
//		Polynomial x("x0");
//		Polynomial y("x1");
//
//		const auto center_node = this->center_node();
//		const auto x_c = center_node.at(0);
//		const auto y_c = center_node.at(1);
//
//		//1 (x - x_c) (y - y_c)  ...
//		for (ushort a = 0; a <= polynomial_order; ++a)
//			for (ushort b = 0; b <= a; ++b)
//				initial_basis_set[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
//	}
//	else if (space_dimension == 3) {
//		Polynomial x("x0");
//		Polynomial y("x1");
//		Polynomial z("x2");
//
//		const auto center_node = this->center_node();
//		const auto x_c = center_node.at(0);
//		const auto y_c = center_node.at(1);
//		const auto z_c = center_node.at(2);
//
//		//1 (x - x_c) (y - y_c) (z - z_c) ...
//		for (ushort a = 0; a <= polynomial_order; ++a)
//			for (ushort b = 0; b <= a; ++b)
//				for (ushort c = 0; c <= b; ++c)
//					initial_basis_set[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ (b - c)) * ((z - z_c) ^ c);
//	}
//	else
//		throw std::runtime_error("not supported space dimension");
//
//	Vector_Function<Polynomial, num_basis> initial_basis_function = initial_basis_set;
//	return initial_basis_function;
//}
//

//template <ushort polynomial_degree>
//auto Geometry::orthonormal_basis_vector_function(void) const {
	//const auto initial_basis_set = this->initial_basis_function<polynomial_degree>();
	//return ms::Gram_Schmidt_process(initial_basis_set, *this);
//}
//
//
////auto Geometry::initial_basis_function(const ushort polynomial_order) const {
	//constexpr auto num_basis = ms::combination_with_repetition(1 + space_dimension, polynomial_order);

	//std::array<Polynomial, num_basis> initial_basis_set = { 0 };
	//
	//ushort index = 0;
	//if (space_dimension == 2) {
	//	Polynomial x("x0");
	//	Polynomial y("x1");

	//	const auto center_node = this->center_node();
	//	const auto x_c = center_node.at(0);
	//	const auto y_c = center_node.at(1);

	//	//1 (x - x_c) (y - y_c)  ...
	//	for (ushort a = 0; a <= polynomial_order; ++a)
	//		for (ushort b = 0; b <= a; ++b)
	//			initial_basis_set[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
	//}
	//else if (space_dimension == 3) {
	//	Polynomial x("x0");
	//	Polynomial y("x1");
	//	Polynomial z("x2");

	//	const auto center_node = this->center_node();
	//	const auto x_c = center_node.at(0);
	//	const auto y_c = center_node.at(1);
	//	const auto z_c = center_node.at(2);

	//	//1 (x - x_c) (y - y_c) (z - z_c) ...
	//	for (ushort a = 0; a <= polynomial_order; ++a)
	//		for (ushort b = 0; b <= a; ++b)
	//			for (ushort c = 0; c <= b; ++c)
	//				initial_basis_set[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ (b - c)) * ((z - z_c) ^ c);
	//}
	//else
	//	throw std::runtime_error("not supported space dimension");

	//Vector_Function<Polynomial, num_basis> initial_basis_function = initial_basis_set;
	//return initial_basis_function;
////}
////
//
////auto Geometry::orthonormal_basis_vector_function(const ushort polynomial_order) const {
////	const auto initial_basis_set = this->initial_basis_function<order>();
////	return ms::Gram_Schmidt_process(initial_basis_set, *this);
////}
//

//const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort integrand_order) const {
//	if (this->integrand_order_to_quadrature_rule_.find(integrand_order) == this->integrand_order_to_quadrature_rule_.end())
//		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->reference_geometry_.quadrature_rule(this->mapping_function_, integrand_order));
//
//	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
//}
//




//
//

//bool Geometry::operator==(const Geometry& other) const {
//	return this->reference_geometry_ == other.reference_geometry_ && this->nodes_ == other.nodes_;
//}
//