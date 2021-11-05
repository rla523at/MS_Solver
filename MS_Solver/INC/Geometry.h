#pragma once
#include "Reference_Geometry.h"

class Geometry
{
public:
	Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes);

public:
	Euclidean_Vector center_node(void) const;
	ushort num_post_nodes(const ushort post_order) const;
	ushort num_post_elements(const ushort post_order) const;
	std::vector<Euclidean_Vector> post_nodes(const ushort post_order) const;
	std::vector<std::vector<int>> post_connectivities(const ushort post_order, const size_t connectivity_start_index) const;
	Vector_Function<Polynomial> orthonormal_basis_vector_function(const ushort solution_order) const;	
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;

private:
	Vector_Function<Polynomial> initial_basis_vector_function(const ushort solution_order) const;
	ushort check_space_dimension(void) const;
	Vector_Function<Polynomial> make_mapping_function(void) const;
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const;
	Irrational_Function make_scale_function(const Vector_Function<Polynomial>& mapping_function) const;

private:
	ushort space_dimension_;
	std::unique_ptr<Reference_Geometry> reference_geometry_;
	std::vector<Euclidean_Vector> nodes_;
	Vector_Function<Polynomial> mapping_function_;
	mutable std::map<size_t, Quadrature_Rule> integrand_order_to_quadrature_rule_;


//public:
//	bool operator==(const Geometry & other) const;
//
//public:
//	Space_Vector_ center_node(void) const;
//	std::vector<Space_Vector_> post_nodes(const ushort post_order) const;
//	std::vector<Space_Vector_> vertex_nodes(void) const;
//	double volume(void) const;
//	std::array<double, space_dimension> projected_volume(void) const;
//	Space_Vector_ normalized_normal_vector(const Space_Vector_ & node) const;
//	std::vector<Geometry> face_geometries(void) const;
//	std::vector<Geometry> sub_simplex_geometries(void) const;
//	bool can_be_periodic_pair(const Geometry & other) const;
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
	//bool is_axis_parallel_node(const Space_Vector_ & node) const;
	//bool is_on_same_axis(const Geometry & other) const;
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
	double integrate(const Polynomial& integrand, const Quadrature_Rule& quadrature_rule);	
	double integrate(const Polynomial& integrand, const Geometry& geometry) {
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand.degree());
		return ms::integrate(integrand, quadrature_rule);
	}
	
	double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry) {
		const auto integrand_degree = f1.degree() + f2.degree();

		const auto quadrature_rule = geometry.get_quadrature_rule(integrand_degree);
		const auto& QP_set = quadrature_rule.nodes;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += f1(QP_set[i]) * f2(QP_set[i]) * QW_set[i];

		return result;
	}
	
	double L2_Norm(const Polynomial& function, const Geometry& geometry);
	Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& functions, const Geometry& geometry);
}