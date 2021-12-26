#pragma once
#include "Euclidean_Vector.h"
#include "Matrix.h"
#include "Polynomial.h"

#include <array>
#include <map>

using uint = unsigned int;

enum class Figure
{
	point = 0,
	line = 1,
	triangle = 2, quadrilateral = 3,
	tetrahedral = 4, hexahedral = 5, prism = 6, pyramid = 7,
	num_figures = 8
};

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

class Reference_Geometry
{
public://Query
	bool operator==(const Reference_Geometry& other) const;
	bool operator!=(const Reference_Geometry& other) const;

	const std::vector<Euclidean_Vector>& get_mapping_nodes(void) const;
	const Vector_Function<Polynomial>& get_mapping_monomial_vector_function(void) const;
	const Matrix& get_inverse_mapping_monomial_matrix(void) const;
	const std::vector<Euclidean_Vector>& get_post_points(const ushort post_order) const;
	const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
	std::vector<ushort> vertex_node_index_sequneces(void) const;

	virtual Euclidean_Vector center_point(void) const abstract;
	virtual std::vector<std::unique_ptr<Reference_Geometry>> face_reference_geometries(void) const abstract;
	virtual Figure figure(void) const abstract;
	virtual bool is_simplex(void) const abstract;
	virtual bool is_line(void) const abstract;
	virtual ushort num_vertices(void) const abstract;
	virtual ushort num_post_nodes(const ushort post_order) const abstract;
	virtual ushort num_post_elements(const ushort post_order) const abstract;
	virtual Vector_Function<Polynomial> normal_vector_function(const Vector_Function<Polynomial>& mapping_function) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_vertex_index_sequences(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const abstract;
	virtual std::vector<std::unique_ptr<Reference_Geometry>> sub_simplex_reference_geometries(void) const abstract;
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
	virtual std::vector<Euclidean_Vector> make_mapping_nodes(void) const abstract;
	virtual std::vector<Euclidean_Vector> make_post_nodes(const ushort post_order) const abstract;
	virtual std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const abstract;
	virtual Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const abstract;

protected:
	ushort order_ = 0;

	static constexpr int num_figures_ = static_cast<int>(Figure::num_figures);

	static inline std::vector<std::vector<std::vector<Euclidean_Vector>>> figure_index_and_order_to_mapping_points_table_ = std::vector<std::vector<std::vector<Euclidean_Vector>>>(num_figures_); //[figure_index][figure_order]
	static inline std::vector<std::vector<Vector_Function<Polynomial>>> figure_index_and_order_to_mapping_monomial_vf_ = std::vector<std::vector<Vector_Function<Polynomial>>>(num_figures_);; //[figure_index][figure_order]
	static inline std::vector<std::vector<Matrix>> figure_index_and_order_to_inverse_mapping_monomial_m_ = std::vector<std::vector<Matrix>>(num_figures_); //[figure_index][figure_order]
	static inline std::vector<std::vector<std::vector<Euclidean_Vector>>> set_of_post_nodes_ = std::vector<std::vector<std::vector<Euclidean_Vector>>>(num_figures_); //figure_index][post_order]
	static inline std::vector<std::vector<std::vector<std::vector<uint>>>> set_of_connectivities_ = std::vector<std::vector<std::vector<std::vector<uint>>>>(num_figures_); //[figure_index][post_order]
	static inline std::vector<std::vector<Quadrature_Rule>> quadrature_rules_ = std::vector<std::vector<Quadrature_Rule>>(num_figures_); //[figure_index][integrand_order]

	//array give incomprehensible error on test get_quadrature_rule_2,4,6
	//static inline std::array<std::vector<std::vector<Euclidean_Vector>>, num_figures_> set_of_mapping_nodes_; //[figure_index][figure_order]
	//static inline std::array<std::vector<Vector_Function<Polynomial>>, num_figures_> set_of_mapping_monomial_vector_function_; //[figure_index][figure_order]
	//static inline std::array<std::vector<Matrix>, num_figures_> set_of_inverse_mapping_monomial_matrix_; //[figure_index][figure_order]
	//static inline std::array<std::vector<std::vector<Euclidean_Vector>>, num_figures_> set_of_post_nodes_; //figure_index][post_order]
	//static inline std::array<std::vector<std::vector<std::vector<uint>>>, num_figures_> set_of_connectivities_; //[figure_index][post_order]
	//static inline std::array<std::vector<Quadrature_Rule>, num_figures_> quadrature_rules_; //[figure_index][integrand_order]
};