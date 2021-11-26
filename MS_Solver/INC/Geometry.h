#pragma once
#include "Reference_Geometry.h"

class Geometry
{
public:
	Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes); // �������� ������ ������ ����� --> �����ʿ���� ������ ?
	Geometry(std::unique_ptr<Reference_Geometry>&& reference_goemetry, std::vector<Euclidean_Vector>&& consisting_nodes);

public://Command
	void change_nodes(std::vector<Euclidean_Vector>&& nodes);

public://Query
	bool operator==(const Geometry& other) const;

	Euclidean_Vector center_point(void) const;
	std::vector<Geometry> face_geometries(void) const;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
	bool is_line(void) const;
	ushort num_post_nodes(const ushort post_order) const;
	ushort num_post_elements(const ushort post_order) const;
	Euclidean_Vector normalized_normal_vector(const Euclidean_Vector& node) const;
	std::vector<Euclidean_Vector> normalized_normal_vectors(const std::vector<Euclidean_Vector>& points) const;
	Vector_Function<Polynomial> orthonormal_basis_vector_function(const ushort solution_order) const;
	std::vector<Euclidean_Vector> post_nodes(const ushort post_order) const;
	std::vector<Euclidean_Vector> post_element_centers(const ushort post_order) const;
	std::vector<std::vector<int>> post_connectivities(const ushort post_order, const size_t connectivity_start_index) const;
	std::vector<double> projected_volume(void) const;
	std::vector<std::vector<Euclidean_Vector>> set_of_face_nodes(void) const;
	double volume(void) const;

protected:	
	ushort check_space_dimension(void) const;
	Vector_Function<Polynomial> initial_basis_vector_function(const ushort solution_order) const;
	bool is_axis_parallel_node(const Euclidean_Vector& node) const;
	bool is_on_axis_plane(const ushort axis_tag) const;
	bool is_on_this_axis_plane(const ushort axis_tag, const double reference_value) const;
	bool is_on_same_axis_plane(const Geometry& other) const;
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
//	std::vector<Geometry> sub_simplex_geometries(void) const;
};


template <typename Function>
Vector_Function<Function> operator*(const Matrix& matrix, const Vector_Function<Function>& vector_function) 
{
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