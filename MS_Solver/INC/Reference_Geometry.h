#pragma once
#include "Euclidean_Vector.h"
#include "Figure.h"
#include "Matrix.h"
#include "Polynomial.h"


#include <array>
#include <map>

using uint = unsigned int;



struct Quadrature_Rule
{
	std::vector<Euclidean_Vector> points;
	std::vector<double> weights;

	bool operator==(const Quadrature_Rule& other) const
	{
		return this->points == other.points && this->weights == other.weights;
	}
	bool is_empty(void) const
	{
		return this->points.empty();
	}
};

class Reference_Geometry abstract
{
public://Query
	bool operator==(const Reference_Geometry& other) const;
	bool operator!=(const Reference_Geometry& other) const;

	ushort check_space_dimension(const std::vector<Euclidean_Vector>& points) const;
	const std::vector<Euclidean_Vector>& get_post_points(const ushort post_order) const;
	const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
	Vector_Function<Polynomial> make_mapping_function(const std::vector<Euclidean_Vector>& points) const;
	std::vector<ushort> vertex_node_index_sequneces(void) const;

	virtual Euclidean_Vector center_point(void) const abstract;
	virtual std::vector<std::shared_ptr<const Reference_Geometry>> face_reference_geometries(void) const abstract;
	virtual Figure figure(void) const abstract;
	virtual bool is_simplex(void) const abstract;
	virtual bool is_line(void) const abstract;
	virtual ushort num_vertices(void) const abstract;
	virtual ushort num_post_nodes(const ushort post_order) const abstract;
	virtual ushort num_post_elements(const ushort post_order) const abstract;
	virtual Vector_Function<Polynomial> make_normal_vector_function(const Vector_Function<Polynomial>& mapping_function) const abstract;
	virtual Euclidean_Vector random_point(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_vertex_index_sequences(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const abstract;
	virtual std::vector<std::shared_ptr<const Reference_Geometry>> sub_simplex_reference_geometries(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_sub_simplex_vertex_index_sequences(void) const abstract;
	virtual Irrational_Function scale_function(const Vector_Function<Polynomial>& mapping_function) const abstract;
	virtual ushort scale_function_order(void) const abstract;

protected:
	void initialize(void);

	Matrix make_inverse_mapping_monomial_matrix(void) const;
	std::vector<std::vector<uint>> quadrilateral_connectivities(const std::array<uint, 4>& node_indexes) const;
	//std::vector<std::vector<uint>> sliced_hexahedral_connectivities(const std::array<uint, 7>& node_indexes) const;
	//std::vector<std::vector<uint>> hexahedral_connectivities(const std::array<uint, 8>& node_indexes) const;

	virtual Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const abstract;
	virtual std::vector<Euclidean_Vector> make_mapping_points(void) const abstract;
	virtual std::vector<Euclidean_Vector> make_post_points(const ushort post_order) const abstract;
	virtual std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const abstract;
	virtual Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const abstract;

protected:
	ushort order_ = 0;

	std::map<ushort, std::vector<Euclidean_Vector>> order_to_mapping_points_;
	std::map<ushort, Vector_Function<Polynomial>> order_to_mapping_monomial_vf_;
	std::map<ushort, Matrix> order_to_inverse_mapping_monomial_m_;
	mutable std::map<ushort, std::vector<Euclidean_Vector>> post_order_to_post_points_;
	mutable std::map<ushort, std::vector<std::vector<uint>>> post_order_to_connectivities_;
	mutable std::map<ushort, Quadrature_Rule> integrand_order_to_quadrature_rule_;
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