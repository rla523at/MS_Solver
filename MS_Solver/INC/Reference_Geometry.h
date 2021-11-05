#pragma once
#include "Euclidean_Vector.h"
#include "Matrix.h"
#include "Polynomial.h"

#include <array>
#include <map>

using uint = unsigned int;

enum class Figure
{
	point,
	line,
	triangle, quadrilateral,
	tetrahedral, hexahedral, prism, pyramid,
	not_in_list
};

struct Quadrature_Rule
{
	std::vector<Euclidean_Vector> nodes;
	std::vector<double> weights;

	bool operator==(const Quadrature_Rule& other) const {
		return this->nodes == other.nodes && this->weights == other.weights;
	}

	bool is_empty(void) const {
		return this->nodes.empty();
	}
};

//template <ushort space_dimension>
//ushort ReferenceGeometry<space_dimension>::scale_function_order(void) const {
//	dynamic_require(this->figure_order_ == 1, "high order mesh is not supported yet");
//
//	switch (this->figure_) {
//	case Figure::line:
//	case Figure::triangle:
//	case Figure::tetrahedral:	return 0;
//	case Figure::quadrilateral: return 1;
//	case Figure::hexahedral:	return 2;
//	default:
//		throw std::runtime_error("not supported figure");
//		return NULL;
//	}
//}

class Reference_Geometry
{
public://Query
	virtual Euclidean_Vector center_node(void) const abstract;
	virtual ushort num_vertex(void) const abstract;
	virtual ushort num_post_nodes(const ushort post_order) const abstract;
	virtual ushort num_post_elements(const ushort post_order) const abstract;
	virtual bool is_simplex(void) const abstract;
	virtual bool is_line(void) const abstract;
	virtual std::vector<ushort> vertex_node_index_sequneces(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_vertex_node_index_sequences(void) const abstract;
	virtual std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const abstract;
	virtual Irrational_Function scale_function(const Vector_Function<Polynomial>& mapping_function) const abstract;
	virtual ushort scale_function_order(void) const abstract;
	virtual const std::vector<Euclidean_Vector>& get_mapping_nodes(void) const abstract;
	virtual const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const abstract;
	virtual const std::vector<Euclidean_Vector>& get_post_nodes(const ushort post_order) const abstract;
	virtual const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const abstract;
	virtual const Vector_Function<Polynomial>& get_mapping_monomial_vector_function(void) const abstract;
	virtual const Matrix& get_inverse_mapping_monomial_matrix(void) const abstract;

protected:
	virtual Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const abstract;
	virtual std::vector<Euclidean_Vector> make_mapping_nodes(void) const abstract;
	virtual std::vector<Euclidean_Vector> make_post_nodes(const ushort post_order) const abstract;
	virtual std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const abstract;
	virtual Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const abstract;

	Matrix make_inverse_mapping_monomial_matrix(void) const;
	std::vector<std::vector<uint>> quadrilateral_connectivities(const std::array<uint, 4>& node_indexes) const;
	std::vector<std::vector<uint>> sliced_hexahedral_connectivities(const std::array<uint, 7>& node_indexes) const;
	std::vector<std::vector<uint>> hexahedral_connectivities(const std::array<uint, 8>& node_indexes) const;

protected:
	ushort order_ = 0;

	//virtual std::vector<ReferenceGeometry> face_reference_geometries(void) const abstract;
	//Vector_Function<Polynomial<space_dimension>, space_dimension> normal_vector_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const;
	//virtual Irrational_Function<space_dimension> scale_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const abstract;
	//virtual std::vector<ReferenceGeometry> sub_simplex_reference_geometries(void) const abstract;
	//virtual std::vector<std::vector<ushort>> set_of_sub_simplex_vertex_node_index_orders(void) const abstract;
	//virtual ushort scale_function_order(void) const abstract;
};

class Reference_Line : public Reference_Geometry
{
public:
	Reference_Line(const ushort order);

public://Query
	Euclidean_Vector center_node(void) const override;
	ushort num_vertex(void) const override;
	ushort num_post_nodes(const ushort post_order) const override;
	ushort num_post_elements(const ushort post_order) const override;
	bool is_simplex(void) const override;
	bool is_line(void) const override;
	std::vector<ushort> vertex_node_index_sequneces(void) const override;
	std::vector<std::vector<ushort>> set_of_face_vertex_node_index_sequences(void) const override;
	std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const override;
	Irrational_Function scale_function(const Vector_Function<Polynomial>& mapping_function) const override;
	ushort scale_function_order(void) const override;
	const std::vector<Euclidean_Vector>& get_mapping_nodes(void) const override;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const override;
	const std::vector<Euclidean_Vector>& get_post_nodes(const ushort post_order) const override;
	const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const override;
	const Vector_Function<Polynomial>& get_mapping_monomial_vector_function(void) const override;
	const Matrix& get_inverse_mapping_monomial_matrix(void) const override;

private:
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const override;
	std::vector<Euclidean_Vector> make_mapping_nodes(void) const override;
	std::vector<Euclidean_Vector> make_post_nodes(const ushort post_order) const override;
	std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const override;
	Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const override;

private:
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_mapping_nodes_;
	static inline std::vector<Quadrature_Rule> quadrature_rules_;
	static inline std::vector<Vector_Function<Polynomial>> set_of_mapping_monomial_vector_function_;
	static inline std::vector<Matrix> set_of_inverse_mapping_monomial_matrix_;
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes_;
	static inline std::vector<std::vector<std::vector<uint>>> set_of_connectivities_;
};

class Reference_Triangle : public Reference_Geometry
{
public:
	Reference_Triangle(const ushort order);

public://Query
	Euclidean_Vector center_node(void) const override;
	ushort num_vertex(void) const override;
	ushort num_post_nodes(const ushort post_order) const override;
	ushort num_post_elements(const ushort post_order) const override;
	bool is_simplex(void) const override;
	bool is_line(void) const override;
	std::vector<ushort> vertex_node_index_sequneces(void) const override;
	std::vector<std::vector<ushort>> set_of_face_vertex_node_index_sequences(void) const override;
	std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const override;
	Irrational_Function scale_function(const Vector_Function<Polynomial>& mapping_function) const override;
	ushort scale_function_order(void) const override;
	const std::vector<Euclidean_Vector>& get_mapping_nodes(void) const override;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const override;
	const std::vector<Euclidean_Vector>& get_post_nodes(const ushort post_order) const override;
	const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const override;
	const Vector_Function<Polynomial>& get_mapping_monomial_vector_function(void) const override;
	const Matrix& get_inverse_mapping_monomial_matrix(void) const override;

private:
	std::vector<Euclidean_Vector> make_mapping_nodes(void) const override;
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const override;
	std::vector<Euclidean_Vector> make_post_nodes(const ushort post_order) const override;
	std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const override;
	Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const override;

private:
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_mapping_nodes_;
	static inline std::vector<Quadrature_Rule> quadrature_rules_;
	static inline std::vector<Vector_Function<Polynomial>> set_of_mapping_monomial_vector_function_;
	static inline std::vector<Matrix> set_of_inverse_mapping_monomial_matrix_;
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes_;
	static inline std::vector<std::vector<std::vector<uint>>> set_of_connectivities_;
};

class Reference_Quadrilateral : public Reference_Geometry
{
public:
	Reference_Quadrilateral(ushort order);

public://Query
	Euclidean_Vector center_node(void) const override;
	ushort num_vertex(void) const override;
	ushort num_post_nodes(const ushort post_order) const override;
	ushort num_post_elements(const ushort post_order) const override;
	bool is_simplex(void) const override;
	bool is_line(void) const override;
	std::vector<ushort> vertex_node_index_sequneces(void) const override;
	std::vector<std::vector<ushort>> set_of_face_vertex_node_index_sequences(void) const override;
	std::vector<std::vector<ushort>> set_of_face_node_index_sequences(void) const override;
	Irrational_Function scale_function(const Vector_Function<Polynomial>& mapping_function) const override;
	ushort scale_function_order(void) const override;
	const std::vector<Euclidean_Vector>& get_mapping_nodes(void) const override;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const override;
	const std::vector<Euclidean_Vector>& get_post_nodes(const ushort post_order) const override;
	const std::vector<std::vector<uint>>& get_connectivities(const ushort post_order) const override;
	const Vector_Function<Polynomial>& get_mapping_monomial_vector_function(void) const override;
	const Matrix& get_inverse_mapping_monomial_matrix(void) const override;

private:
	std::vector<Euclidean_Vector> make_mapping_nodes(void) const override;
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const override;
	std::vector<Euclidean_Vector> make_post_nodes(const ushort post_order) const override;
	std::vector<std::vector<uint>> make_connectivities(const ushort post_order) const override;
	Vector_Function<Polynomial> make_mapping_monomial_vector_function(void) const override;

private:
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_mapping_nodes_;
	static inline std::vector<Quadrature_Rule> quadrature_rules_;
	static inline std::vector<Vector_Function<Polynomial>> set_of_mapping_monomial_vector_function_;
	static inline std::vector<Matrix> set_of_inverse_mapping_monomial_matrix_;
	static inline std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes_;
	static inline std::vector<std::vector<std::vector<uint>>> set_of_connectivities_;
};

class Reference_Geometry_Factory
{
public:
	static std::unique_ptr<Reference_Geometry> make(const Figure figure, const ushort order) {
		switch (figure)
		{
		case Figure::line:			return std::make_unique<Reference_Line>(order);
		case Figure::triangle:		return std::make_unique<Reference_Triangle>(order);
		case Figure::quadrilateral: return std::make_unique<Reference_Quadrilateral>(order);
		default: EXCEPTION("not supproted figure"); return nullptr;
		}
	}
};
