#pragma once
#include "Matrix.h"
#include "Polynomial.h"

#include <map>
#include <set>
#include <unordered_set>


using uint		= unsigned int;


enum class Figure
{
	point,
	line,
	triangle, quadrilateral,
	tetrahedral, hexahedral, prism, pyramid,
	not_in_list
};


enum class ElementType
{
	cell, inner_face,
	slip_wall_2D,
	supersonic_outlet_2D,
	x_periodic, y_periodic,
	not_in_list
};


enum class FaceType
{
	inward_face,
	outward_face,
	not_my_face
};


template <ushort space_dimension>
struct Quadrature_Rule
{
	std::vector<Euclidean_Vector<space_dimension>> points;
	std::vector<double> weights;

	bool operator==(const Quadrature_Rule& other) const {
		return this->points == other.points && this->weights == other.weights;
	}
};


template <ushort space_dimension>
class ReferenceGeometry
{
private:
	using Space_Vector_ = Euclidean_Vector<space_dimension>;

private:
	Figure figure_;
	ushort figure_order_;

	inline static std::map<std::pair<Figure, ushort>, std::vector<Space_Vector_>> key_to_mapping_nodes_;
	inline static std::map<std::pair<Figure, ushort>, Dynamic_Vector_Function_<Polynomial<space_dimension>>> key_to_mapping_monomial_vector_function_;
	inline static std::map<std::pair<Figure, ushort>, Dynamic_Matrix_> key_to_inverse_mapping_monomial_matrix_;
	inline static std::map<std::pair<Figure, ushort>, Quadrature_Rule<space_dimension>> key_to_reference_quadrature_rule_;
	inline static std::map<std::pair<Figure, ushort>, std::vector<Space_Vector_>> key_to_reference_post_nodes_;
	inline static std::map<std::pair<Figure, ushort>, std::vector<std::vector<ushort>>> key_to_reference_connectivity_;

public:
	ReferenceGeometry(const Figure figure, const ushort figure_order);

	bool operator==(const ReferenceGeometry& other) const;
	bool operator!=(const ReferenceGeometry& other) const;

	Space_Vector_ center_node(void) const;
	Vector_Function<Polynomial<space_dimension>, space_dimension> normal_vector_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const;
	ushort num_vertex(void) const;
	std::vector<ushort> vertex_node_index_orders(void) const;
	std::vector<std::vector<ushort>> face_vertex_node_index_orders_set(void) const;
	std::vector<std::vector<ushort>> face_node_index_orders_set(void) const;
	std::vector<ReferenceGeometry> face_reference_geometries(void) const;
	//std::vector<std::vector<ushort>> local_connectivities(void) const;
	Vector_Function<Polynomial<space_dimension>, space_dimension> mapping_function(const std::vector<Space_Vector_>& mapped_nodes) const;
	Irrational_Function<space_dimension> scale_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const;
	Quadrature_Rule<space_dimension> quadrature_rule(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function, const ushort physical_integrand_order) const;
	std::vector<Space_Vector_> post_nodes(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function, const ushort post_order) const;
	std::vector<std::vector<size_t>> post_connectivities(const ushort post_order, const size_t connectivity_start_index) const;

	//private: for test
	std::vector<Space_Vector_> mapping_nodes(void) const;
	Dynamic_Vector_Function_<Polynomial<space_dimension>> mapping_monomial_vector_function(void) const;
	Dynamic_Matrix_ inverse_mapping_monomial_matrix(void) const;
	Quadrature_Rule<space_dimension> reference_quadrature_rule(const ushort integrand_order) const;
	std::vector<Space_Vector_> reference_post_nodes(const ushort post_order) const;
	std::vector<std::vector<ushort>> reference_connectivity(const ushort post_order) const;
};


//template <ushort space_dimension>
//Dynamic_Vector_Function_<Polynomial<space_dimension>> operator*(const Dynamic_Matrix_& A, const Dynamic_Vector_Function_<Polynomial<space_dimension>>& v);


template <ushort space_dimension>
class Geometry
{
private:
	using Space_Vector_ = Euclidean_Vector<space_dimension>;

public:
	ReferenceGeometry<space_dimension> reference_geometry_;

private:
	std::vector<Space_Vector_> nodes_;
	Vector_Function<Polynomial<space_dimension>, space_dimension> mapping_function_;
	mutable std::map<size_t, Quadrature_Rule<space_dimension>> integrand_order_to_quadrature_rule_;

public:
	Geometry(const ReferenceGeometry<space_dimension> reference_geometry, std::vector<Space_Vector_>&& consisting_nodes)
		: reference_geometry_(reference_geometry), nodes_(std::move(consisting_nodes)), mapping_function_(this->reference_geometry_.mapping_function(this->nodes_)) {};
		
	Space_Vector_ center_node(void) const;
	Space_Vector_ normalized_normal_vector(const Space_Vector_& node) const;
	double volume(void) const;
	std::array<double, space_dimension> coordinate_projected_volume(void) const;
	std::vector<Geometry> faces_geometry(void) const;
	std::vector<Space_Vector_> post_nodes(const ushort post_order) const;
	std::vector<Space_Vector_> vertex_nodes(void) const;
	bool is_axis_parallel(const Geometry& other, const ushort axis_tag) const;
	const Quadrature_Rule<space_dimension>& get_quadrature_rule(const ushort integrand_order) const;

	template <ushort order>
	auto initial_basis_function(void) const;

	template <ushort order>
	auto orthonormal_basis_function(void) const;

	//private: for test
	std::vector<std::vector<Space_Vector_>> calculate_faces_nodes(void) const;
	bool is_axis_parallel_node(const Space_Vector_& node, const size_t axis_tag) const;
};


template <ushort space_dimension>
class Element
{
public:
	Geometry<space_dimension> geometry_;

private:
	ElementType element_type_;
	std::vector<uint> node_indexes_;

public:
	Element(const ElementType element_type, Geometry<space_dimension>&& geometry, std::vector<uint>&& node_indexes)
		: element_type_(element_type), geometry_(std::move(geometry)), node_indexes_(std::move(node_indexes)) {};

	ElementType type(void) const;
	Euclidean_Vector<space_dimension> normalized_normal_vector(const Element& owner_cell_element, const Euclidean_Vector<space_dimension>& node) const;
	std::vector<uint> vertex_node_indexes(void) const;
	std::vector<Element> make_inner_face_elements(void) const;
	bool is_periodic_pair(const Element& other) const;
	std::vector<std::pair<uint, uint>> find_periodic_vnode_index_pairs(const Element& other) const;
	std::vector<std::vector<uint>> face_node_indexes_set(void) const;
	std::vector<std::vector<uint>> face_vertex_node_indexes_set(void) const;

	//private:
	bool is_periodic_boundary(void) const;
	FaceType check_face_type(const Element& owner_cell_element) const;
};


namespace ms {
	template <ushort space_dimension>
	double integrate(const Polynomial<space_dimension>& integrand, const Quadrature_Rule<space_dimension>& quadrature_rule) {
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += integrand(QP_set[i]) * QW_set[i];

		return result;
	}

	template <ushort space_dimension>
	double integrate(const Polynomial<space_dimension>& integrand, const Geometry<space_dimension>& geometry) {
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand.order());
		return ms::integrate(integrand, quadrature_rule);
	}

	template <ushort space_dimension>
	double inner_product(const Polynomial<space_dimension>& f1, const Polynomial<space_dimension>& f2, const Geometry<space_dimension>& geometry) {
		return ms::integrate(f1 * f2, geometry);
	}

	template <ushort space_dimension>
	double L2_Norm(const Polynomial<space_dimension>& function, const Geometry<space_dimension>& geometry) {
		return std::sqrt(ms::inner_product(function, function, geometry));
	}

	template <ushort space_dimension, ushort range_dimension>
	Vector_Function<Polynomial<space_dimension>, range_dimension> Gram_Schmidt_process(const Vector_Function<Polynomial<space_dimension>, range_dimension>& function, const Geometry<space_dimension>& geometry) {				
		std::array<Polynomial<space_dimension>, range_dimension> normalized_function = { 0 };

		for (ushort i = 0; i < range_dimension; ++i) {
			normalized_function[i] = function[i];

			for (ushort j = 0; j < i; ++j)
				normalized_function[i] -= ms::inner_product(normalized_function[i], normalized_function[j], geometry) * normalized_function[j];

			normalized_function[i] *= 1.0 / ms::L2_Norm(normalized_function[i], geometry);
		}

		return normalized_function;
	}

	template <ushort space_dimension>
	Dynamic_Vector_Function_<Polynomial<space_dimension>> Gram_Schmidt_process(const Dynamic_Vector_Function_<Polynomial<space_dimension>>& functions, const Geometry<space_dimension>& geometry) {
		const auto range_dimension = functions.range_dimension();

		std::vector<Polynomial<space_dimension>> normalized_functions(range_dimension);

		for (ushort i = 0; i < range_dimension; ++i) {
			normalized_functions[i] = functions[i];

			for (ushort j = 0; j < i; ++j)
				normalized_functions[i] -= ms::inner_product(normalized_functions[i], normalized_functions[j], geometry) * normalized_functions[j];

			normalized_functions[i] *= 1.0 / ms::L2_Norm(normalized_functions[i], geometry);
		}

		return normalized_functions;
	}

	template <ushort space_dimension>
	void gemv(const Dynamic_Matrix_& A, const Dynamic_Vector_Function_<Polynomial<space_dimension>>& v, Polynomial<space_dimension>* ptr) {
		const auto [num_row, num_column] = A.size();
		const auto range_dimension = v.range_dimension();
		dynamic_require(num_column == range_dimension, "number of column should be same with range dimension");

		for (size_t i = 0; i < num_row; ++i)
			for (size_t j = 0; j < num_column; ++j)
				ptr[i] += A.at(i, j) * v.at(j);
	}

	//double inner_product(const Polynomial& f1, const Polynomial& f2, const QuadratureRule& quadrature_rule);
	//double L2_Norm(const Polynomial& polynomial, const QuadratureRule& quadrature_rule);
	//std::vector<Polynomial> Gram_Schmidt_Process(const std::vector<Polynomial>& initial_polynomial_set, const QuadratureRule& quadrature_rule);
}


//template definition part
template <ushort space_dimension>
ReferenceGeometry<space_dimension>::ReferenceGeometry(const Figure figure, const ushort figure_order) : figure_(figure), figure_order_(figure_order) {
	dynamic_require(figure_order > 0, "figure order should be greater than 0");

	auto key = std::make_pair(this->figure_, this->figure_order_);
	if (ReferenceGeometry::key_to_mapping_nodes_.find(key) == ReferenceGeometry::key_to_mapping_nodes_.end()) {
		ReferenceGeometry::key_to_mapping_nodes_.emplace(key, this->mapping_nodes());
		ReferenceGeometry::key_to_mapping_monomial_vector_function_.emplace(key, this->mapping_monomial_vector_function());
		ReferenceGeometry::key_to_inverse_mapping_monomial_matrix_.emplace(std::move(key), this->inverse_mapping_monomial_matrix());
	}
};


template <ushort space_dimension>
bool ReferenceGeometry<space_dimension>::operator==(const ReferenceGeometry& other) const {
	return this->figure_ == other.figure_ && this->figure_order_ == other.figure_order_;
}

template <ushort space_dimension>
bool ReferenceGeometry<space_dimension>::operator != (const ReferenceGeometry& other) const {
	return !((*this) == other);
}

template <ushort space_dimension>
Euclidean_Vector<space_dimension> ReferenceGeometry<space_dimension>::center_node(void) const {
	if constexpr (space_dimension == 2) {
		switch (this->figure_) {
		case Figure::line:			return { 0, 0 };
		case Figure::triangle:		return { -1.0 / 3.0, -1.0 / 3.0 };
		case Figure::quadrilateral:	return { 0, 0 };
		default:
			throw std::runtime_error("wrong figure");
			return {};
		}
	}
	else if constexpr (space_dimension == 3) {
		switch (this->figure_) {
		case Figure::triangle:		return { -1.0 / 3.0, -1.0 / 3.0, 0 };
		case Figure::quadrilateral:	return { 0, 0, 0 };
		default:
			throw std::runtime_error("wrong figure");
			return {};
		}
	}
	else {
		throw std::runtime_error("unsupported space dimension");
		return {};
	}
}

template <ushort space_dimension>
Vector_Function<Polynomial<space_dimension>, space_dimension> ReferenceGeometry<space_dimension>::normal_vector_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const {
	if constexpr (space_dimension == 2) {
		dynamic_require(this->figure_ == Figure::line, "In 2D case, normal vector must be determined only for line");

		constexpr ushort r = 0;
		const auto T_r = mapping_function.differentiate<r>();

		std::array<Polynomial<space_dimension>, space_dimension> normal_vector_function = { 0 };
		normal_vector_function[0] = -1 * T_r.at(1);
		normal_vector_function[1] = T_r.at(0);

		return normal_vector_function;
	}
	else if constexpr (space_dimension == 3) {
		dynamic_require(this->figure_ == Figure::triangle || this->figure_ == Figure::quadrilateral, "In 3D case, normal vector must be determined only for triangle or quadrilateral");

		constexpr ushort r = 0;
		constexpr ushort s = 1;
		const auto T_r = mapping_function.differentiate<r>();
		const auto T_s = mapping_function.differentiate<s>();

		return T_r.cross_product(T_s);
	}
	else {
		throw std::runtime_error("not supported figure");
		return {};
	}
}


template <ushort space_dimension>
ushort ReferenceGeometry<space_dimension>::num_vertex(void) const {
	switch (this->figure_) {
	case Figure::line:			return 2;
	case Figure::triangle:		return 3;
	case Figure::quadrilateral:	return 4;
	default:
		throw std::runtime_error("wrong element figure");
		return NULL;
	}
}

template <ushort space_dimension>
std::vector<ushort> ReferenceGeometry<space_dimension>::vertex_node_index_orders(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		return { 0,1 };
	}
	case Figure::triangle: {
		//  2
		//  弛 \
		//	弛  \
		//  0式式式1
		return { 0,1,2 };
	}
	case Figure::quadrilateral: {
		//  3式式式式式2
		//  弛     弛
		//  0式式式式式1
		return { 0,1,2,3 };
	}
	default:
		throw std::runtime_error("wrong element figure");
		return { 0 };
	}
}

template <ushort space_dimension>
std::vector<std::vector<ushort>> ReferenceGeometry<space_dimension>::face_vertex_node_index_orders_set(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const std::vector<ushort> face0_node_index = { 0 };
		const std::vector<ushort> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		const std::vector<ushort> face0_node_index = { 0,1 };
		const std::vector<ushort> face1_node_index = { 1,2 };
		const std::vector<ushort> face2_node_index = { 2,0 };
		return { face0_node_index,face1_node_index, face2_node_index };
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		std::vector<ushort> face0_node_index = { 0,1 };
		std::vector<ushort> face1_node_index = { 1,2 };
		std::vector<ushort> face2_node_index = { 2,3 };
		std::vector<ushort> face3_node_index = { 3,0 };
		return { face0_node_index,face1_node_index, face2_node_index,face3_node_index };
	}
	default:
		throw std::runtime_error("wrong element figure");
		return std::vector<std::vector<ushort>>();
	}
}

template <ushort space_dimension>
std::vector<std::vector<ushort>> ReferenceGeometry<space_dimension>::face_node_index_orders_set(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const std::vector<ushort> face0_node_index = { 0 };
		const std::vector<ushort> face1_node_index = { 1 };
		return { face0_node_index,face1_node_index };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		constexpr ushort num_face = 3;
		std::vector<std::vector<ushort>> face_node_index_orders_set(num_face);
		face_node_index_orders_set[0] = { 0,1 };
		face_node_index_orders_set[1] = { 1,2 };
		face_node_index_orders_set[2] = { 2,0 };

		if (this->figure_order_ > 1) {
			const auto num_additional_point = this->figure_order_ - 1;

			ushort index = num_face;
			for (ushort iface = 0; iface < num_face; ++iface)
				for (ushort ipoint = 0; ipoint < num_additional_point; ++ipoint)
					face_node_index_orders_set[iface].push_back(index++);
		}

		return face_node_index_orders_set;
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		constexpr ushort num_face = 4;
		std::vector<std::vector<ushort>> face_node_index_orders_set(num_face);
		face_node_index_orders_set[0] = { 0,1 };
		face_node_index_orders_set[1] = { 1,2 };
		face_node_index_orders_set[2] = { 2,3 };
		face_node_index_orders_set[3] = { 3,0 };

		if (this->figure_order_ > 1) {
			const ushort num_additional_point = this->figure_order_ - 1;

			ushort index = num_face;
			for (ushort iface = 0; iface < num_face; ++iface)
				for (ushort ipoint = 0; ipoint < num_additional_point; ++ipoint)
					face_node_index_orders_set[iface].push_back(index++);
		}

		return face_node_index_orders_set;

	}
	default:
		throw std::runtime_error("wrong element figure");
		return std::vector<std::vector<ushort>>();
	}
}

template <ushort space_dimension>
std::vector<ReferenceGeometry<space_dimension>> ReferenceGeometry<space_dimension>::face_reference_geometries(void) const {
	switch (this->figure_) {
	case Figure::line: {
		// 0 式式式式 1
		const ReferenceGeometry face0_reference_geometry = { Figure::point,this->figure_order_ };
		const ReferenceGeometry face1_reference_geometry = { Figure::point,this->figure_order_ };
		return { face0_reference_geometry,face1_reference_geometry };
	}
	case Figure::triangle: {
		//      2
		//  2  / \  1
		//	  /   \
		//   0式式式式式1
		//      0
		const ReferenceGeometry face0_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face1_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face2_refrence_geometry = { Figure::line,this->figure_order_ };
		return { face0_refrence_geometry,face1_refrence_geometry, face2_refrence_geometry };
	}
	case Figure::quadrilateral: {
		//      2
		//   3式式式式式2
		//3  弛     弛   1
		//   0式式式式式1
		//      0
		const ReferenceGeometry face0_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face1_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face2_refrence_geometry = { Figure::line,this->figure_order_ };
		const ReferenceGeometry face3_refrence_geometry = { Figure::line,this->figure_order_ };
		return { face0_refrence_geometry,face1_refrence_geometry, face2_refrence_geometry,face3_refrence_geometry };
	}
	default:
		throw std::runtime_error("not supported figure");
		return std::vector<ReferenceGeometry>();
	}
}

template <ushort space_dimension>
Vector_Function<Polynomial<space_dimension>, space_dimension> ReferenceGeometry<space_dimension>::mapping_function(const std::vector<Space_Vector_>& mapped_nodes) const {
	const auto key = std::make_pair(this->figure_, this->figure_order_);
	const auto& mapping_nodes = ReferenceGeometry::key_to_mapping_nodes_.at(key);

	const auto num_mapping_node = mapping_nodes.size();
	const auto num_mapped_node = mapped_nodes.size();
	dynamic_require(num_mapping_node == num_mapped_node, "number of mapping node should be same with mapped node");

	//	X = CM
	//	X : mapped node matrix			
	//	C : mapping coefficient matrix	
	//	M : mapping monomial matrix
	Dynamic_Matrix_ X(space_dimension, num_mapped_node);
	for (size_t j = 0; j < num_mapped_node; ++j)
		X.change_column(j, mapped_nodes[j]);

	const auto& inv_M = ReferenceGeometry::key_to_inverse_mapping_monomial_matrix_.at(key);
	const auto C = X * inv_M;

	const auto& monomial_vector_function = ReferenceGeometry::key_to_mapping_monomial_vector_function_.at(key);

	Vector_Function<Polynomial<space_dimension>, space_dimension> mapping_function;
	ms::gemv(C, monomial_vector_function, mapping_function.data());
	return mapping_function;
}

template <ushort space_dimension>
Irrational_Function<space_dimension> ReferenceGeometry<space_dimension>::scale_function(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function) const {
	switch (this->figure_)
	{
	case Figure::line: {
		constexpr size_t r = 0;
		const auto T_r = mapping_function.differentiate<r>();
		return T_r.L2_norm();
	}
	case Figure::triangle:
	case Figure::quadrilateral: {
		constexpr size_t r = 0;
		constexpr size_t s = 1;
		const auto T_r = mapping_function.differentiate<r>();
		const auto T_s = mapping_function.differentiate<s>();
		const auto cross_product = T_r.cross_product(T_s);
		return cross_product.L2_norm();
	}
	default:
		throw std::runtime_error("not supported figure");
		return Irrational_Function<space_dimension>();
	}
}

template <ushort space_dimension>
Quadrature_Rule<space_dimension> ReferenceGeometry<space_dimension>::quadrature_rule(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function, const ushort physical_integrand_order) const {
	const auto scale_function = this->scale_function(mapping_function);

	constexpr ushort heuristic_additional_order = 2;
	const auto reference_integrand_order = physical_integrand_order + heuristic_additional_order;
	const auto key = std::make_pair(this->figure_, reference_integrand_order);
	if (ReferenceGeometry::key_to_reference_quadrature_rule_.find(key) == ReferenceGeometry::key_to_reference_quadrature_rule_.end())
		ReferenceGeometry::key_to_reference_quadrature_rule_.emplace(key, this->reference_quadrature_rule(reference_integrand_order));

	const auto& reference_quadrature_rule = ReferenceGeometry::key_to_reference_quadrature_rule_.at(key);

	const auto transformed_QP = mapping_function(reference_quadrature_rule.points);

	const auto num_QP = transformed_QP.size();
	std::vector<double> transformed_QW(num_QP);

	for (size_t i = 0; i < num_QP; ++i) {
		const auto& point = reference_quadrature_rule.points[i];
		const auto& weight = reference_quadrature_rule.weights[i];

		transformed_QW[i] = scale_function(point) * weight;
	}

	return { transformed_QP, transformed_QW };
}


template <ushort space_dimension>
std::vector<Euclidean_Vector<space_dimension>> ReferenceGeometry<space_dimension>::post_nodes(const Vector_Function<Polynomial<space_dimension>, space_dimension>& mapping_function, const ushort post_order) const {
	const auto key = std::make_pair(this->figure_, post_order);
	if (ReferenceGeometry::key_to_reference_post_nodes_.find(key) == ReferenceGeometry::key_to_reference_post_nodes_.end())
		ReferenceGeometry::key_to_reference_post_nodes_.emplace(key, this->reference_post_nodes(post_order));

	const auto& reference_post_nodes = ReferenceGeometry::key_to_reference_post_nodes_.at(key);
	return mapping_function(reference_post_nodes);
}

template <ushort space_dimension>
std::vector<std::vector<size_t>> ReferenceGeometry<space_dimension>::post_connectivities(const ushort post_order, const size_t connectivity_start_index) const {
	const auto key = std::make_pair(this->figure_, post_order);
	if (ReferenceGeometry::key_to_reference_connectivity_.find(key) == ReferenceGeometry::key_to_reference_connectivity_.end())
		ReferenceGeometry::key_to_reference_connectivity_.emplace(key, this->reference_connectivity(post_order));

	const auto& reference_connectivities = ReferenceGeometry::key_to_reference_connectivity_.at(key);

	const auto num_connectivity = reference_connectivities.size();
	
	std::vector<std::vector<size_t>> connectivities(num_connectivity);	

	for (ushort i = 0; i < num_connectivity; ++i) {
		auto& connectivity = connectivities[i];
		const auto& reference_connecitivity = reference_connectivities[i];

		const auto num_index = reference_connecitivity.size();
		connectivity.resize(num_index);

		for (ushort j = 0; j < num_index; ++j)
			connectivity[j] = reference_connecitivity[j] + connectivity_start_index;
	}

	return connectivities;
}

//template <ushort space_dimension>
//std::vector<std::vector<ushort>> ReferenceGeometry<space_dimension>::local_connectivities(void) const {
//	switch (this->figure_) {
//	case Figure::triangle: {
//		//  2
//		//  弛 \ 
//		//  弛  \
//		//  0式式式1
//		return { { 0,1,2 } };
//	}
//	case Figure::quadrilateral: {
//		//  3式式式式式2
//		//  弛     弛
//		//  弛     弛 
//		//  0式式式式式1
//
//		return { {0,1,2},{0,2,3} };
//	}
//	default:
//		throw std::runtime_error("wrong element figure");
//		return {};
//	}
//}

template <ushort space_dimension>
std::vector<Euclidean_Vector<space_dimension>> ReferenceGeometry<space_dimension>::mapping_nodes(void) const {
	if constexpr (space_dimension == 2) {
		switch (this->figure_) {
		case Figure::line: {
			switch (this->figure_order_)
			{
			case 1: {
				// 0 式式式式 1 		
				return { { -1,0 }, { 1,0 } };
			}
			case 2: {
				// 0 式式 2 式式 1
				return { { -1,0 }, { 1,0 }, { 0,0 } };
			}
			default:	throw std::runtime_error("unsuported figure order");
			}
			break;
		}
		case Figure::triangle: {
			switch (this->figure_order_)
			{
			case 1: {
				//	  2
				//    弛 \ 
				//    弛  \  
				//	  弛   \
				//    0式式式式式1
				return { { -1, -1 }, { 1, -1 }, { -1, 1 } };
			}
			default:	throw std::runtime_error("unsuported figure order");
			}
			break;
		}
		case Figure::quadrilateral: {
			switch (this->figure_order_)
			{
			case 1: {
				//   3式式式式式2
				//   弛     弛   
				//   0式式式式式1
				return { { -1, -1 }, { 1, -1 }, { 1, 1 }, { -1, 1 } };
			}
			default:	throw std::runtime_error("unsuported figure order");
			}
			break;
		}
		default:
			throw std::runtime_error("not supported element figure");
			return {};
		}
	}
}

template <ushort space_dimension>
Dynamic_Vector_Function_<Polynomial<space_dimension>> ReferenceGeometry<space_dimension>::mapping_monomial_vector_function(void) const {
	Polynomial<space_dimension> r("x0");
	Polynomial<space_dimension> s("x1");

	switch (this->figure_) {
	case Figure::line: {
		const auto num_monomial = this->figure_order_ + 1;
		std::vector<Polynomial<space_dimension>> mapping_monomial_vector(num_monomial);

		for (ushort a = 0, index = 0; a <= this->figure_order_; ++a)
			mapping_monomial_vector[index++] = (r ^ a);

		return mapping_monomial_vector;	// 1 r r^2 ...
		//return std::move(mapping_monomial_vector);	// 1 r r^2 ...
	}
	case Figure::triangle: {
		const auto num_monomial = static_cast<size_t>((this->figure_order_ + 2) * (this->figure_order_ + 1) * 0.5);
		std::vector<Polynomial<space_dimension>> mapping_monomial_vector(num_monomial);

		for (ushort a = 0, index = 0; a <= this->figure_order_; ++a)
			for (ushort b = 0; b <= a; ++b)
				mapping_monomial_vector[index++] = (r ^ (a - b)) * (s ^ b);

		return mapping_monomial_vector;	// 1 r s r^2 rs s^2 ...
	}
	case Figure::quadrilateral: {
		const auto num_monomial = (this->figure_order_ + 1) * (this->figure_order_ + 1);
		std::vector<Polynomial<space_dimension>> mapping_monomial_vector(num_monomial);

		for (ushort a = 0, index = 0; a <= this->figure_order_; ++a) {
			for (ushort b = 0; b <= a; ++b)
				mapping_monomial_vector[index++] = (r ^ a) * (s ^ b);

			if (a == 0)
				continue;

			for (int c = static_cast<int>(a - 1); 0 <= c; --c)
				mapping_monomial_vector[index++] = (r ^ c) * (s ^ a);
		}
		return mapping_monomial_vector;	// 1 r rs s r^2 r^2s r^2s^2 rs^2 s^2...
	}
	default:
		throw std::runtime_error("not supproted figure");
		return std::vector<Polynomial<space_dimension>>(NULL);
	}
}

template <ushort space_dimension>
Dynamic_Matrix_ ReferenceGeometry<space_dimension>::inverse_mapping_monomial_matrix(void) const {
	const auto key = std::make_pair(this->figure_, this->figure_order_);
	const auto mapping_nodes = ReferenceGeometry::key_to_mapping_nodes_.at(key);
	const auto mapping_monomial_vector_function = ReferenceGeometry::key_to_mapping_monomial_vector_function_.at(key);

	const auto matrix_order = mapping_monomial_vector_function.range_dimension();
	Dynamic_Matrix_ transformation_monomial_matrix(matrix_order);
	for (size_t i = 0; i < matrix_order; ++i)
		transformation_monomial_matrix.change_column(i, mapping_monomial_vector_function(mapping_nodes[i]));

	return transformation_monomial_matrix.be_inverse();
}

template <ushort space_dimension>
Quadrature_Rule<space_dimension> ReferenceGeometry<space_dimension>::reference_quadrature_rule(const ushort integrand_order) const {
	switch (this->figure_)
	{
	case Figure::line: {
		switch (integrand_order)
		{
		case 0:
		case 1:		return { { { 0.000000000000000 } }, { 2.000000000000000 } };
		case 2:
		case 3:		return { { { -0.577350269189626 }, { 0.577350269189626 } }, { 1.000000000000000, 1.000000000000000 } };
		case 4:
		case 5:		return { { { -0.774596669241483 }, { 0.000000000000000 }, { 0.774596669241483 } }, { 0.555555555555554, 0.888888888888889, 0.555555555555554 } };
		case 6:
		case 7:		return { { { -0.861136311594052 }, { -0.339981043584856 }, { 0.339981043584856 }, { 0.861136311594052 } }, { 0.347854845137454, 0.652145154862546, 0.652145154862546, 0.347854845137454 } };
		case 8:
		case 9:		return { { { -0.906179845938664 }, { -0.538469310105683 }, { 0.000000000000000 }, { 0.538469310105683 }, { 0.906179845938664 } }, { 0.236926885056189, 0.478628670499366, 0.568888888888889, 0.478628670499366, 0.236926885056189 } };
		case 10:
		case 11:	return { { { -0.932469514203152 }, { -0.661209386466264 }, { -0.238619186083197 }, { 0.238619186083197 }, { 0.661209386466264 }, { 0.932469514203152 } }, { 0.171324492379171, 0.360761573048139, 0.467913934572691, 0.467913934572691, 0.360761573048139, 0.171324492379171 } };
		case 12:
		case 13:	return { { { -0.949107912342758 }, { -0.741531185599394 }, { -0.405845151377397 }, { 0.000000000000000 }, { 0.405845151377397 }, { 0.741531185599394 }, { 0.949107912342758 } }, { 0.129484966168870, 0.279705391489277, 0.381830050505119, 0.417959183673469, 0.381830050505119, 0.279705391489277, 0.129484966168870 } };
		case 14:
		case 15:	return { { { -0.960289856497536 }, { -0.796666477413627 }, { -0.525532409916329 }, { -0.183434642495650 }, { 0.183434642495650 }, { 0.525532409916329 }, { 0.796666477413627 }, { 0.960289856497536 } }, { 0.101228536290377, 0.222381034453374, 0.313706645877887, 0.362683783378362, 0.362683783378362, 0.313706645877887, 0.222381034453374, 0.101228536290377 } };
		case 16:
		case 17:	return { { { -0.968160239507626 }, { -0.836031107326636 }, { -0.613371432700590 }, { -0.324253423403809 }, { 0.000000000000000 }, { 0.324253423403809 }, { 0.613371432700590 }, { 0.836031107326636 }, { 0.968160239507626 } }, { 0.081274388361575, 0.180648160694857, 0.260610696402936, 0.312347077040003, 0.330239355001260, 0.312347077040003, 0.260610696402936, 0.180648160694857, 0.081274388361575 } };
		case 18:
		case 19:	return { { { -0.973906528517172 }, { -0.865063366688985 }, { -0.679409568299024 }, { -0.433395394129247 }, { -0.148874338981631 }, { 0.148874338981631 }, { 0.433395394129247 }, { 0.679409568299024 }, { 0.865063366688985 }, { 0.973906528517172 } }, { 0.066671344308688, 0.149451349150581, 0.219086362515982, 0.269266719309996, 0.295524224714753, 0.295524224714753, 0.269266719309996, 0.219086362515982, 0.149451349150581, 0.066671344308688 } };
		case 20:
		case 21:	return { { { -0.978228658146057 }, { -0.887062599768095 }, { -0.730152005574049 }, { -0.519096129206812 }, { -0.269543155952345 }, { 0.000000000000000 }, { 0.269543155952345 }, { 0.519096129206812 }, { 0.730152005574049 }, { 0.887062599768095 }, { 0.978228658146057 } }, { 0.055668567116174, 0.125580369464904, 0.186290210927734, 0.233193764591990, 0.262804544510247, 0.272925086777901, 0.262804544510247, 0.233193764591990, 0.186290210927734, 0.125580369464904, 0.055668567116174 } };
		default:
			throw std::runtime_error("not supported order");
		}
	}
	case Figure::triangle: {
		switch (integrand_order)
		{
		case 0:		return { { { -0.5, 0 } }, { 2 } };
		case 1:
		case 2:		return { { { -0.6666666666666667, -0.577350269189626 }, { -0.9106836025229592, 0.577350269189626 }, { 0.2440169358562927, -0.577350269189626 }, { -0.6666666666666667, 0.577350269189626 } }, { 0.788675134594813, 0.211324865405187, 0.788675134594813, 0.211324865405187 } };
		case 3:
		case 4:		return { { { -0.7999999999999996, -0.774596669241483 }, { -0.8872983346207415, 0 }, { -0.9745966692414834, 0.774596669241483 }, { -0.1127016653792585, -0.774596669241483 }, { -0.5, 0 }, { -0.8872983346207415, 0.774596669241483 }, { 0.5745966692414826, -0.774596669241483 }, { -0.1127016653792585, 0 }, { -0.7999999999999996, 0.774596669241483 } }, { 0.2738575106854125, 0.2469135802469129, 0.03478446462322775, 0.4381720170966613, 0.3950617283950618, 0.05565514339716456, 0.2738575106854125, 0.2469135802469129, 0.03478446462322775 } };
		case 5:
		case 6:		return { { { -0.8707778735729041, -0.861136311594052 }, { -0.9069626449468777, -0.339981043584856 }, { -0.9541736666471743, 0.339981043584856 }, { -0.9903584380211479, 0.861136311594052 }, { -0.3858073769376817, -0.861136311594052 }, { -0.5577935549985239, -0.339981043584856 }, { -0.7821874885863321, 0.339981043584856 }, { -0.9541736666471743, 0.861136311594052 }, { 0.2469436885317338, -0.861136311594052 }, { -0.1022254014166201, -0.339981043584856 }, { -0.5577935549985238, 0.339981043584856 }, { -0.9069626449468777, 0.861136311594052 }, { 0.7319141851669562, -0.861136311594052 }, { 0.2469436885317338, -0.339981043584856 }, { -0.3858073769376817, 0.339981043584856 }, { -0.8707778735729041, 0.861136311594052 } }, { 0.1126015323077027, 0.1519885905918008, 0.07486326126005106, 0.008401460977899435, 0.211101109416918, 0.2849424819989601, 0.140350821011734, 0.01575074243493392, 0.211101109416918, 0.2849424819989601, 0.140350821011734, 0.01575074243493392, 0.1126015323077027, 0.1519885905918008, 0.07486326126005106, 0.008401460977899435 } };
		case 7:
		case 8:		return { { { -0.9105809565927103, -0.906179845938664 }, { -0.9278302861536236, -0.538469310105683 }, { -0.9530899229693319, 0 }, { -0.9783495597850402, 0.538469310105683 }, { -0.9955988893459535, 0.906179845938664 }, { -0.5601197503206428, -0.906179845938664 }, { -0.644974598962845, -0.538469310105683 }, { -0.7692346550528415, 0 }, { -0.893494711142838, 0.538469310105683 }, { -0.9783495597850402, 0.906179845938664 }, { -0.04691007703066802, -0.906179845938664 }, { -0.2307653449471585, -0.538469310105683 }, { -0.5, 0 }, { -0.7692346550528415, 0.538469310105683 }, { -0.9530899229693319, 0.906179845938664 }, { 0.4662995962593067, -0.906179845938664 }, { 0.1834439090685281, -0.538469310105683 }, { -0.2307653449471585, 0 }, { -0.644974598962845, 0.538469310105683 }, { -0.9278302861536237, 0.906179845938664 }, { 0.8167608025313744, -0.906179845938664 }, { 0.4662995962593067, -0.538469310105683 }, { -0.04691007703066802, 0 }, { -0.5601197503206428, 0.538469310105683 }, { -0.9105809565927103, 0.906179845938664 } }, { 0.05350108223322567, 0.08723120988299211, 0.06739253619376044, 0.02616879011700774, 0.002633266629202917, 0.1080803972647221, 0.1762204318958822, 0.1361432662753753, 0.05286497232810846, 0.005319602735277746, 0.1284622942592381, 0.2094522369422109, 0.1618172839506173, 0.06283429560853965, 0.006322778128282769, 0.1080803972647221, 0.1762204318958822, 0.1361432662753753, 0.05286497232810846, 0.005319602735277746, 0.05350108223322567, 0.08723120988299211, 0.06739253619376044, 0.02616879011700774, 0.002633266629202917 } };
		case 9:
		case 10:	return { { { -0.9347496974591312, -0.9324695142031521 }, { -0.9439088615608247, -0.661209386466264 }, { -0.9581777223232526, -0.238619186083197 }, { -0.9742917918798994, 0.238619186083197 }, { -0.9885606526423274, 0.661209386466264 }, { -0.9977198167440209, 0.9324695142031521 }, { -0.6726487338239366, -0.9324695142031521 }, { -0.7185989263755467, -0.661209386466264 }, { -0.7901837230061085, -0.238619186083197 }, { -0.8710256634601555, 0.238619186083197 }, { -0.9426104600907174, 0.661209386466264 }, { -0.9885606526423274, 0.9324695142031521 }, { -0.2643273942032976, -0.9324695142031521 }, { -0.3675935226230415, -0.661209386466264 }, { -0.5284695579835037, -0.238619186083197 }, { -0.7101496280996933, 0.238619186083197 }, { -0.8710256634601555, 0.661209386466264 }, { -0.9742917918798994, 0.9324695142031521 }, { 0.1967969084064496, -0.9324695142031521 }, { 0.0288029090893055, -0.661209386466264 }, { -0.2329112559332993, -0.238619186083197 }, { -0.5284695579835037, 0.238619186083197 }, { -0.7901837230061085, 0.661209386466264 }, { -0.9581777223232526, 0.9324695142031521 }, { 0.6051182480270888, -0.9324695142031521 }, { 0.3798083128418107, -0.661209386466264 }, { 0.0288029090893055, -0.238619186083197 }, { -0.3675935226230415, 0.238619186083197 }, { -0.7185989263755467, 0.661209386466264 }, { -0.9439088615608248, 0.9324695142031521 }, { 0.8672192116622832, -0.9324695142031521 }, { 0.6051182480270888, -0.661209386466264 }, { 0.1967969084064496, -0.238619186083197 }, { -0.2643273942032976, 0.238619186083197 }, { -0.6726487338239368, 0.661209386466264 }, { -0.9347496974591312, 0.9324695142031521 } }, { 0.02836100152117781, 0.0513374279511389, 0.049647026182223, 0.03051809113558391, 0.01046986542124473, 0.0009910801678028134, 0.05972035509877095, 0.1081022976149208, 0.1045427831942518, 0.06426258389333622, 0.02204661497324695, 0.002086938273612684, 0.0774583226595905, 0.1402105301458922, 0.1355937790222319, 0.08334967114506463, 0.02859483694169573, 0.002706794658216405, 0.0774583226595905, 0.1402105301458922, 0.1355937790222319, 0.08334967114506463, 0.02859483694169573, 0.002706794658216405, 0.05972035509877095, 0.1081022976149208, 0.1045427831942518, 0.06426258389333622, 0.02204661497324695, 0.002086938273612684, 0.02836100152117781, 0.0513374279511389, 0.049647026182223, 0.03051809113558391, 0.01046986542124473, 0.0009910801678028134 } };
		case 11:
		case 12:	return { { { -0.9504029146358142, -0.949107912342758 }, { -0.9556849211223275, -0.741531185599394 }, { -0.9642268026617964, -0.405845151377397 }, { -0.974553956171379, 0 }, { -0.9848811096809615, 0.405845151377397 }, { -0.9934229912204304, 0.741531185599394 }, { -0.9987049977069438, 0.949107912342758 }, { -0.7481081943789636, -0.949107912342758 }, { -0.7749342496082214, -0.741531185599394 }, { -0.8183164352463219, -0.405845151377397 }, { -0.870765592799697, 0 }, { -0.9232147503530721, 0.405845151377397 }, { -0.9665969359911726, 0.741531185599394 }, { -0.9934229912204304, 0.949107912342758 }, { -0.4209640416964355, -0.949107912342758 }, { -0.4826304010243249, -0.741531185599394 }, { -0.5823551434482712, -0.405845151377397 }, { -0.7029225756886985, 0 }, { -0.8234900079291259, 0.405845151377397 }, { -0.9232147503530722, 0.741531185599394 }, { -0.9848811096809615, 0.949107912342758 }, { -0.02544604382862098, -0.949107912342758 }, { -0.129234407200303, -0.741531185599394 }, { -0.2970774243113015, -0.405845151377397 }, { -0.5, 0 }, { -0.7029225756886985, 0.405845151377397 }, { -0.870765592799697, 0.741531185599394 }, { -0.974553956171379, 0.949107912342758 }, { 0.3700719540391935, -0.949107912342758 }, { 0.2241615866237189, -0.741531185599394 }, { -0.01179970517433182, -0.405845151377397 }, { -0.2970774243113015, 0 }, { -0.5823551434482711, 0.405845151377397 }, { -0.8183164352463218, 0.741531185599394 }, { -0.9642268026617964, 0.949107912342758 }, { 0.6972161067217216, -0.949107912342758 }, { 0.5164654352076156, -0.741531185599394 }, { 0.2241615866237189, -0.405845151377397 }, { -0.129234407200303, 0 }, { -0.4826304010243249, 0.405845151377397 }, { -0.7749342496082215, 0.741531185599394 }, { -0.9556849211223276, 0.949107912342758 }, { 0.8995108269785723, -0.949107912342758 }, { 0.6972161067217215, -0.741531185599394 }, { 0.3700719540391935, -0.405845151377397 }, { -0.02544604382862098, 0 }, { -0.4209640416964355, 0.405845151377397 }, { -0.7481081943789636, 0.741531185599394 }, { -0.9504029146358142, 0.949107912342758 } }, { 0.01633971902233046, 0.03153707751100931, 0.03475337161903316, 0.02705971537896383, 0.01468787955288011, 0.004680565643230263, 0.0004266374414229519, 0.03529604741916743, 0.06812443847836634, 0.07507207749196593, 0.05845271854796313, 0.03172784626693878, 0.01011066754980336, 0.0009215957350721348, 0.04818316692765089, 0.09299769892288511, 0.1024820257759689, 0.07979468810555948, 0.04331216169277281, 0.0138022248360196, 0.001258084244262363, 0.05274230535088141, 0.1017972322343419, 0.1121789753788724, 0.08734493960849631, 0.04741040083224651, 0.01510820486158434, 0.001377125407046246, 0.04818316692765089, 0.09299769892288511, 0.1024820257759689, 0.07979468810555948, 0.04331216169277281, 0.0138022248360196, 0.001258084244262363, 0.03529604741916743, 0.06812443847836634, 0.07507207749196593, 0.05845271854796313, 0.03172784626693878, 0.01011066754980336, 0.0009215957350721348, 0.01633971902233046, 0.03153707751100931, 0.03475337161903316, 0.02705971537896383, 0.01468787955288011, 0.004680565643230263, 0.0004266374414229519 } };
		case 13:
		case 14:	return { { { -0.9610783042460291, -0.960289856497536 }, { -0.9643270581779191, -0.796666477413627 }, { -0.9697104445422814, -0.525532409916329 }, { -0.9765028202603553, -0.18343464249565 }, { -0.9837870362371807, 0.18343464249565 }, { -0.9905794119552546, 0.525532409916329 }, { -0.9959627983196169, 0.796666477413627 }, { -0.9992115522515068, 0.960289856497536 }, { -0.8007036790940101, -0.960289856497536 }, { -0.8173387381173185, -0.796666477413627 }, { -0.844904060636017, -0.525532409916329 }, { -0.8796840326953073, -0.18343464249565 }, { -0.9169824447183197, 0.18343464249565 }, { -0.95176241677761, 0.525532409916329 }, { -0.9793277392963085, 0.796666477413627 }, { -0.9959627983196169, 0.960289856497536 }, { -0.5349529979610744, -0.960289856497536 }, { -0.573769993138719, -0.796666477413627 }, { -0.6380921569362322, -0.525532409916329 }, { -0.719249308576779, -0.18343464249565 }, { -0.8062831013395499, 0.18343464249565 }, { -0.8874402529800968, 0.525532409916329 }, { -0.95176241677761, 0.796666477413627 }, { -0.9905794119552546, 0.960289856497536 }, { -0.1996476062584693, -0.960289856497536 }, { -0.2664521977773303, -0.796666477413627 }, { -0.3771515411561002, -0.525532409916329 }, { -0.5168241340337535, -0.18343464249565 }, { -0.6666105084618966, 0.18343464249565 }, { -0.8062831013395499, 0.525532409916329 }, { -0.9169824447183198, 0.796666477413627 }, { -0.9837870362371808, 0.960289856497536 }, { 0.1599374627560052, -0.960289856497536 }, { 0.06311867519095721, -0.796666477413627 }, { -0.0973160489275709, -0.525532409916329 }, { -0.2997412234705965, -0.18343464249565 }, { -0.5168241340337535, 0.18343464249565 }, { -0.719249308576779, 0.525532409916329 }, { -0.8796840326953073, 0.796666477413627 }, { -0.9765028202603552, 0.960289856497536 }, { 0.4952428544586104, -0.960289856497536 }, { 0.370436470552346, -0.796666477413627 }, { 0.1636245668525612, -0.525532409916329 }, { -0.0973160489275709, -0.18343464249565 }, { -0.3771515411561001, 0.18343464249565 }, { -0.6380921569362322, 0.525532409916329 }, { -0.844904060636017, 0.796666477413627 }, { -0.9697104445422814, 0.960289856497536 }, { 0.760993535591546, -0.960289856497536 }, { 0.6140052155309454, -0.796666477413627 }, { 0.370436470552346, -0.525532409916329 }, { 0.06311867519095721, -0.18343464249565 }, { -0.2664521977773303, 0.18343464249565 }, { -0.573769993138719, 0.525532409916329 }, { -0.8173387381173185, 0.796666477413627 }, { -0.9643270581779191, 0.960289856497536 }, { 0.9213681607435651, -0.960289856497536 }, { 0.760993535591546, -0.796666477413627 }, { 0.4952428544586104, -0.525532409916329 }, { 0.1599374627560052, -0.18343464249565 }, { -0.1996476062584693, 0.18343464249565 }, { -0.5349529979610744, 0.525532409916329 }, { -0.8007036790940101, 0.796666477413627 }, { -0.9610783042460291, 0.960289856497536 } }, { 0.01004375733945304, 0.02022265498028209, 0.02422245286926617, 0.02172427927521026, 0.01498966925243749, 0.007533611717515962, 0.002288651636172856, 0.00020345922003913, 0.02206434300837125, 0.0444255651490272, 0.05321240752324866, 0.04772436582622505, 0.0329296291009185, 0.01655000090197412, 0.005027759335525518, 0.0004469636080836971, 0.03112554564587481, 0.06266989030061787, 0.07506524072180071, 0.06732341526693157, 0.04645289793099652, 0.02334661894615327, 0.00709251812460491, 0.0006305189409073175, 0.03598499044536026, 0.07245416447754377, 0.08678472663211513, 0.07783421639230392, 0.05370531033333868, 0.02699158656581295, 0.008199830449599781, 0.0007289580822874853, 0.03598499044536026, 0.07245416447754377, 0.08678472663211513, 0.07783421639230392, 0.05370531033333868, 0.02699158656581295, 0.008199830449599781, 0.0007289580822874853, 0.03112554564587481, 0.06266989030061787, 0.07506524072180071, 0.06732341526693157, 0.04645289793099652, 0.02334661894615327, 0.00709251812460491, 0.0006305189409073175, 0.02206434300837125, 0.0444255651490272, 0.05321240752324866, 0.04772436582622505, 0.0329296291009185, 0.01655000090197412, 0.005027759335525518, 0.0004469636080836971, 0.01004375733945304, 0.02022265498028209, 0.02422245286926617, 0.02172427927521026, 0.01498966925243749, 0.007533611717515962, 0.002288651636172856, 0.00020345922003913 } };
		case 15:
		case 16:	return { { { -0.9686671246817318, -0.968160239507626 }, { -0.9707706046430857, -0.836031107326636 }, { -0.9743153199987874, -0.61337143270059 }, { -0.9789180440838081, -0.324253423403809 }, { -0.9840801197538129, 0 }, { -0.9892421954238177, 0.324253423403809 }, { -0.9938449195088385, 0.61337143270059 }, { -0.9973896348645401, 0.836031107326636 }, { -0.9994931148258941, 0.968160239507626 }, { -0.8386414724620959, -0.968160239507626 }, { -0.8494740062089006, -0.836031107326636 }, { -0.8677286363546227, -0.61337143270059 }, { -0.891431816272783, -0.324253423403809 }, { -0.918015553663318, 0 }, { -0.944599291053853, 0.324253423403809 }, { -0.9683024709720133, 0.61337143270059 }, { -0.9865571011177354, 0.836031107326636 }, { -0.9973896348645401, 0.968160239507626 }, { -0.6195265131917514, -0.968160239507626 }, { -0.6450689617285768, -0.836031107326636 }, { -0.6881122572265872, -0.61337143270059 }, { -0.7440028980840231, -0.324253423403809 }, { -0.806685716350295, 0 }, { -0.8693685346165668, 0.324253423403809 }, { -0.9252591754740027, 0.61337143270059 }, { -0.9683024709720132, 0.836031107326636 }, { -0.9938449195088385, 0.968160239507626 }, { -0.3350112279799912, -0.968160239507626 }, { -0.379654132349956, -0.836031107326636 }, { -0.4548848887872422, -0.61337143270059 }, { -0.552570141294545, -0.324253423403809 }, { -0.6621267117019045, 0 }, { -0.7716832821092641, 0.324253423403809 }, { -0.8693685346165668, 0.61337143270059 }, { -0.944599291053853, 0.836031107326636 }, { -0.9892421954238177, 0.968160239507626 }, { -0.01591988024618701, -0.968160239507626 }, { -0.081984446336682, -0.836031107326636 }, { -0.193314283649705, -0.61337143270059 }, { -0.3378732882980955, -0.324253423403809 }, { -0.5, 0 }, { -0.6621267117019045, 0.324253423403809 }, { -0.806685716350295, 0.61337143270059 }, { -0.918015553663318, 0.836031107326636 }, { -0.9840801197538129, 0.968160239507626 }, { 0.3031714674876171, -0.968160239507626 }, { 0.215685239676592, -0.836031107326636 }, { 0.06825632148783223, -0.61337143270059 }, { -0.123176435301646, -0.324253423403809 }, { -0.3378732882980955, 0 }, { -0.552570141294545, 0.324253423403809 }, { -0.7440028980840232, 0.61337143270059 }, { -0.891431816272783, 0.836031107326636 }, { -0.9789180440838081, 0.968160239507626 }, { 0.5876867526993774, -0.968160239507626 }, { 0.4811000690552128, -0.836031107326636 }, { 0.3014836899271773, -0.61337143270059 }, { 0.06825632148783223, -0.324253423403809 }, { -0.193314283649705, 0 }, { -0.4548848887872422, 0.324253423403809 }, { -0.6881122572265872, 0.61337143270059 }, { -0.8677286363546228, 0.836031107326636 }, { -0.9743153199987875, 0.968160239507626 }, { 0.8068017119697218, -0.968160239507626 }, { 0.6855051135355366, -0.836031107326636 }, { 0.4811000690552127, -0.61337143270059 }, { 0.215685239676592, -0.324253423403809 }, { -0.081984446336682, 0 }, { -0.379654132349956, 0.324253423403809 }, { -0.6450689617285768, 0.61337143270059 }, { -0.8494740062089006, 0.836031107326636 }, { -0.9707706046430858, 0.968160239507626 }, { 0.9368273641893579, -0.968160239507626 }, { 0.8068017119697217, -0.836031107326636 }, { 0.5876867526993774, -0.61337143270059 }, { 0.3031714674876172, -0.324253423403809 }, { -0.01591988024618701, 0 }, { -0.3350112279799912, 0.324253423403809 }, { -0.6195265131917516, 0.61337143270059 }, { -0.8386414724620959, 0.836031107326636 }, { -0.9686671246817318, 0.968160239507626 } }, { 0.006500367017424581, 0.01347836749000478, 0.01708638995104882, 0.01680862795979199, 0.01342000079532422, 0.008577189683159995, 0.004094584999583912, 0.001203701279113231, 0.0001051591861235364, 0.01444833199254737, 0.02995829738399937, 0.03797783016022493, 0.0373604500255599, 0.02982856603501678, 0.01906447494013144, 0.009101012802371234, 0.002675460578435511, 0.0002337367765706412, 0.02084377636592117, 0.04321911008813613, 0.05478842811273874, 0.05389776935251935, 0.04303195414326739, 0.02750321991429733, 0.01312950696688454, 0.003859732874460045, 0.00033719858471156, 0.02498167846612465, 0.05179895873279031, 0.06566501533832468, 0.06459754318835401, 0.05157464862910972, 0.03296315334707955, 0.01573597392849199, 0.004625966232901031, 0.000404139176827337, 0.02641271197951784, 0.05476617512723753, 0.0694265255080294, 0.06829790500794712, 0.05452901579582412, 0.03485139225027233, 0.01663738277850538, 0.004890956942796017, 0.0004272896111305921, 0.02498167846612465, 0.05179895873279031, 0.06566501533832468, 0.06459754318835401, 0.05157464862910972, 0.03296315334707955, 0.01573597392849199, 0.004625966232901031, 0.000404139176827337, 0.02084377636592117, 0.04321911008813613, 0.05478842811273874, 0.05389776935251935, 0.04303195414326739, 0.02750321991429733, 0.01312950696688454, 0.003859732874460045, 0.00033719858471156, 0.01444833199254737, 0.02995829738399937, 0.03797783016022493, 0.0373604500255599, 0.02982856603501678, 0.01906447494013144, 0.009101012802371234, 0.002675460578435511, 0.0002337367765706412, 0.006500367017424581, 0.01347836749000478, 0.01708638995104882, 0.01680862795979199, 0.01342000079532422, 0.008577189683159995, 0.004094584999583912, 0.001203701279113231, 0.0001051591861235364 } };
		case 17:
		case 18:	return { { { -0.9742469631441846, -0.973906528517172 }, { -0.9756670111138168, -0.865063366688985 }, { -0.9780891871608004, -0.679409568299024 }, { -0.9812988690798357, -0.433395394129247 }, { -0.985010940099215, -0.148874338981631 }, { -0.988895588417957, 0.148874338981631 }, { -0.9926076594373363, 0.433395394129247 }, { -0.9958173413563716, 0.679409568299024 }, { -0.9982395174033551, 0.865063366688985 }, { -0.9996595653729874, 0.973906528517172 }, { -0.8668238492856299, -0.973906528517172 }, { -0.8741673141936407, -0.865063366688985 }, { -0.8866930634517123, -0.679409568299024 }, { -0.903291225656342, -0.433395394129247 }, { -0.9224873823002004, -0.148874338981631 }, { -0.9425759843887846, 0.148874338981631 }, { -0.9617721410326431, 0.433395394129247 }, { -0.9783703032372728, 0.679409568299024 }, { -0.9908960524953444, 0.865063366688985 }, { -0.9982395174033551, 0.973906528517172 }, { -0.6835922269426524, -0.973906528517172 }, { -0.7010392650617513, -0.865063366688985 }, { -0.730798680748133, -0.679409568299024 }, { -0.770233575898957, -0.433395394129247 }, { -0.8158409398478528, -0.148874338981631 }, { -0.8635686284511712, 0.148874338981631 }, { -0.909175992400067, 0.433395394129247 }, { -0.948610887550891, 0.679409568299024 }, { -0.9783703032372727, 0.865063366688985 }, { -0.9958173413563716, 0.973906528517172 }, { -0.4407877346919107, -0.973906528517172 }, { -0.471623253096604, -0.865063366688985 }, { -0.52421940172918, -0.679409568299024 }, { -0.5939157838262227, -0.433395394129247 }, { -0.6745212539831456, -0.148874338981631 }, { -0.7588741401461014, 0.148874338981631 }, { -0.8394796103030243, 0.433395394129247 }, { -0.909175992400067, 0.679409568299024 }, { -0.961772141032643, 0.865063366688985 }, { -0.9926076594373363, 0.973906528517172 }, { -0.1599787505636739, -0.973906528517172 }, { -0.2062983545928464, -0.865063366688985 }, { -0.2853057105304597, -0.679409568299024 }, { -0.3900001988355295, -0.433395394129247 }, { -0.5110817844036087, -0.148874338981631 }, { -0.6377925545780222, 0.148874338981631 }, { -0.7588741401461014, 0.433395394129247 }, { -0.8635686284511712, 0.679409568299024 }, { -0.9425759843887845, 0.865063366688985 }, { -0.988895588417957, 0.973906528517172 }, { 0.1338852790808459, -0.973906528517172 }, { 0.07136172128183138, -0.865063366688985 }, { -0.03528472117051629, -0.679409568299024 }, { -0.1766044070352235, -0.433395394129247 }, { -0.3400438766147603, -0.148874338981631 }, { -0.5110817844036089, 0.148874338981631 }, { -0.6745212539831456, 0.433395394129247 }, { -0.8158409398478528, 0.679409568299024 }, { -0.9224873823002004, 0.865063366688985 }, { -0.985010940099215, 0.973906528517172 }, { 0.4146942632090826, -0.973906528517172 }, { 0.336686619785589, -0.865063366688985 }, { 0.203628970028204, -0.679409568299024 }, { 0.02731117795546967, -0.433395394129247 }, { -0.1766044070352235, -0.148874338981631 }, { -0.3900001988355296, 0.148874338981631 }, { -0.5939157838262227, 0.433395394129247 }, { -0.770233575898957, 0.679409568299024 }, { -0.903291225656342, 0.865063366688985 }, { -0.9812988690798357, 0.973906528517172 }, { 0.6574987554598244, -0.973906528517172 }, { 0.5661026317507363, -0.865063366688985 }, { 0.410208249047157, -0.679409568299024 }, { 0.203628970028204, -0.433395394129247 }, { -0.03528472117051629, -0.148874338981631 }, { -0.2853057105304597, 0.148874338981631 }, { -0.52421940172918, 0.433395394129247 }, { -0.730798680748133, 0.679409568299024 }, { -0.8866930634517123, 0.865063366688985 }, { -0.9780891871608004, 0.973906528517172 }, { 0.8407303778028019, -0.973906528517172 }, { 0.7392306808826257, -0.865063366688985 }, { 0.5661026317507363, -0.679409568299024 }, { 0.3366866197855889, -0.433395394129247 }, { 0.07136172128183144, -0.148874338981631 }, { -0.2062983545928465, 0.148874338981631 }, { -0.471623253096604, 0.433395394129247 }, { -0.7010392650617512, 0.679409568299024 }, { -0.8741673141936406, 0.865063366688985 }, { -0.9756670111138168, 0.973906528517172 }, { 0.9481534916613564, -0.973906528517172 }, { 0.8407303778028018, -0.865063366688985 }, { 0.6574987554598245, -0.679409568299024 }, { 0.4146942632090827, -0.433395394129247 }, { 0.133885279080846, -0.148874338981631 }, { -0.159978750563674, 0.148874338981631 }, { -0.4407877346919107, 0.433395394129247 }, { -0.6835922269426525, 0.679409568299024 }, { -0.8668238492856298, 0.865063366688985 }, { -0.9742469631441845, 0.973906528517172 } }, { 0.004387074522396848, 0.00929185979426592, 0.01226538498559636, 0.01286642521300538, 0.01131813402104741, 0.008384863316467971, 0.005085948940982216, 0.00234139732304471, 0.0006722625623504124, 5.799362953077529e-05, 0.009834123085334445, 0.02082875329379134, 0.02749424588563135, 0.0288415454460565, 0.02537087585156437, 0.01879561823873495, 0.01140072903617321, 0.005248506572875442, 0.001506952469137529, 0.0001299992712818888, 0.01441621147982787, 0.03053365406746336, 0.04030485074533407, 0.04227990792340959, 0.03719212262556008, 0.02755320480255082, 0.01671275815682937, 0.007693983495150223, 0.002209098391043433, 0.0001905708288132014, 0.01771815427246953, 0.03752719596452479, 0.04953642393731129, 0.05196385557058451, 0.0457107449708518, 0.03386409349471978, 0.02054071055738368, 0.009456242142927666, 0.002715078517704924, 0.0002342198815180672, 0.01944593753793903, 0.0411866550814514, 0.05436696119271134, 0.05703110347256456, 0.05016822169208674, 0.03716634570116908, 0.02254373499300701, 0.01037836623539956, 0.002979839008847915, 0.0002570597995763472, 0.01944593753793903, 0.0411866550814514, 0.05436696119271134, 0.05703110347256456, 0.05016822169208674, 0.03716634570116908, 0.02254373499300701, 0.01037836623539956, 0.002979839008847915, 0.0002570597995763472, 0.01771815427246953, 0.03752719596452479, 0.04953642393731129, 0.05196385557058451, 0.0457107449708518, 0.03386409349471978, 0.02054071055738368, 0.009456242142927666, 0.002715078517704924, 0.0002342198815180672, 0.01441621147982787, 0.03053365406746336, 0.04030485074533407, 0.04227990792340959, 0.03719212262556008, 0.02755320480255082, 0.01671275815682937, 0.007693983495150223, 0.002209098391043433, 0.0001905708288132014, 0.009834123085334445, 0.02082875329379134, 0.02749424588563135, 0.0288415454460565, 0.02537087585156437, 0.01879561823873495, 0.01140072903617321, 0.005248506572875442, 0.001506952469137529, 0.0001299992712818888, 0.004387074522396848, 0.00929185979426592, 0.01226538498559636, 0.01286642521300538, 0.01131813402104741, 0.008384863316467971, 0.005085948940982216, 0.00234139732304471, 0.0006722625623504124, 5.799362953077529e-05 } };
		case 19:
		case 20:	return { { { -0.9784656538091175, -0.978228658146057 }, { -0.9794580575203291, -0.887062599768095 }, { -0.9811661346136811, -0.730152005574049 }, { -0.9834636194310185, -0.519096129206812 }, { -0.9861801709767138, -0.269543155952345 }, { -0.9891143290730284, 0 }, { -0.992048487169343, 0.269543155952345 }, { -0.9947650387150384, 0.519096129206812 }, { -0.9970625235323758, 0.730152005574049 }, { -0.9987706006257278, 0.887062599768095 }, { -0.9997630043369393, 0.978228658146057 }, { -0.8882919991423672, -0.978228658146057 }, { -0.8934400279536657, -0.887062599768095 }, { -0.9023005652422252, -0.730152005574049 }, { -0.9142186162325163, -0.519096129206812 }, { -0.9283105482422671, -0.269543155952345 }, { -0.9435312998840475, 0 }, { -0.9587520515258279, 0.269543155952345 }, { -0.9728439835355787, 0.519096129206812 }, { -0.9847620345258697, 0.730152005574049 }, { -0.9936225718144293, 0.887062599768095 }, { -0.9987706006257278, 0.978228658146057 }, { -0.7330894820416732, -0.978228658146057 }, { -0.7453899710481793, -0.887062599768095 }, { -0.7665609756219032, -0.730152005574049 }, { -0.7950374780966583, -0.519096129206812 }, { -0.8287081627645337, -0.269543155952345 }, { -0.8650760027870246, 0 }, { -0.9014438428095154, 0.269543155952345 }, { -0.9351145274773909, 0.519096129206812 }, { -0.963591029952146, 0.730152005574049 }, { -0.9847620345258699, 0.887062599768095 }, { -0.9970625235323759, 0.978228658146057 }, { -0.5243310904917735, -0.978228658146057 }, { -0.5462521456712334, -0.887062599768095 }, { -0.5839816017294213, -0.730152005574049 }, { -0.6347303956787477, -0.519096129206812 }, { -0.6947358910817587, -0.269543155952345 }, { -0.7595480646034061, 0 }, { -0.8243602381250534, 0.269543155952345 }, { -0.8843657335280645, 0.519096129206812 }, { -0.9351145274773909, 0.730152005574049 }, { -0.9728439835355788, 0.887062599768095 }, { -0.9947650387150386, 0.978228658146057 }, { -0.2774946687830019, -0.978228658146057 }, { -0.3107911044265171, -0.887062599768095 }, { -0.3680993131428296, -0.730152005574049 }, { -0.4451829178272916, -0.519096129206812 }, { -0.5363267564603751, -0.269543155952345 }, { -0.6347715779761725, 0 }, { -0.7332163994919698, 0.269543155952345 }, { -0.8243602381250532, 0.519096129206812 }, { -0.9014438428095153, 0.730152005574049 }, { -0.9587520515258279, 0.887062599768095 }, { -0.992048487169343, 0.978228658146057 }, { -0.01088567092697151, -0.978228658146057 }, { -0.05646870011595251, -0.887062599768095 }, { -0.1349239972129755, -0.730152005574049 }, { -0.240451935396594, -0.519096129206812 }, { -0.3652284220238275, -0.269543155952345 }, { -0.5, 0 }, { -0.6347715779761725, 0.269543155952345 }, { -0.7595480646034061, 0.519096129206812 }, { -0.8650760027870246, 0.730152005574049 }, { -0.9435312998840475, 0.887062599768095 }, { -0.9891143290730284, 0.978228658146057 }, { 0.2557233269290589, -0.978228658146057 }, { 0.1978537041946121, -0.887062599768095 }, { 0.0982513187168787, -0.730152005574049 }, { -0.03572095296589628, -0.519096129206812 }, { -0.1941300875872799, -0.269543155952345 }, { -0.3652284220238275, 0 }, { -0.5363267564603751, 0.269543155952345 }, { -0.6947358910817587, 0.519096129206812 }, { -0.8287081627645336, 0.730152005574049 }, { -0.9283105482422671, 0.887062599768095 }, { -0.986180170976714, 0.978228658146057 }, { 0.5025597486378306, -0.978228658146057 }, { 0.4333147454393283, -0.887062599768095 }, { 0.3141336073034702, -0.730152005574049 }, { 0.1538265248855596, -0.519096129206812 }, { -0.03572095296589628, -0.269543155952345 }, { -0.240451935396594, 0 }, { -0.4451829178272917, 0.269543155952345 }, { -0.6347303956787476, 0.519096129206812 }, { -0.7950374780966583, 0.730152005574049 }, { -0.9142186162325163, 0.887062599768095 }, { -0.9834636194310185, 0.978228658146057 }, { 0.7113181401877302, -0.978228658146057 }, { 0.6324525708162743, -0.887062599768095 }, { 0.496712981195952, -0.730152005574049 }, { 0.3141336073034702, -0.519096129206812 }, { 0.0982513187168787, -0.269543155952345 }, { -0.1349239972129755, 0 }, { -0.3680993131428297, 0.269543155952345 }, { -0.5839816017294213, 0.519096129206812 }, { -0.766560975621903, 0.730152005574049 }, { -0.9023005652422251, 0.887062599768095 }, { -0.9811661346136811, 0.978228658146057 }, { 0.8665206572884241, -0.978228658146057 }, { 0.7805026277217607, -0.887062599768095 }, { 0.6324525708162743, -0.730152005574049 }, { 0.4333147454393284, -0.519096129206812 }, { 0.1978537041946121, -0.269543155952345 }, { -0.05646870011595251, 0 }, { -0.3107911044265171, 0.269543155952345 }, { -0.5462521456712334, 0.519096129206812 }, { -0.7453899710481793, 0.730152005574049 }, { -0.8934400279536657, 0.887062599768095 }, { -0.9794580575203291, 0.978228658146057 }, { 0.9566943119551745, -0.978228658146057 }, { 0.8665206572884241, -0.887062599768095 }, { 0.7113181401877302, -0.730152005574049 }, { 0.5025597486378304, -0.519096129206812 }, { 0.2557233269290589, -0.269543155952345 }, { -0.01088567092697151, 0 }, { -0.2774946687830019, 0.269543155952345 }, { -0.5243310904917735, 0.519096129206812 }, { -0.7330894820416731, 0.730152005574049 }, { -0.8882919991423672, 0.887062599768095 }, { -0.9784656538091177, 0.978228658146057 } }, { 0.003065254786336921, 0.006596113363469354, 0.008971278567846241, 0.009860120851096311, 0.009286677986218874, 0.007596674255491598, 0.005343274438285346, 0.003121441884166165, 0.001399230542270533, 0.0003947658625615832, 3.37345784310486e-05, 0.006914778815286163, 0.01487989355803276, 0.02023792843044774, 0.02224303019808677, 0.02094942465785366, 0.0171370166169049, 0.01205366713879896, 0.007041528916287183, 0.003156465085552, 0.0008905356369090305, 7.610041074477408e-05, 0.01025761916059888, 0.02207334252415016, 0.03002163452865245, 0.03299607100161919, 0.03107709234297891, 0.02542163599166264, 0.01788082168660207, 0.01044564459125499, 0.004682408158847183, 0.001321050991849573, 0.0001128899495178914, 0.01284024971520857, 0.02763089812771649, 0.03758038567929434, 0.04130371625698049, 0.03890158328739785, 0.03182221421866715, 0.02238279779882984, 0.01307561558760397, 0.005861329913579828, 0.001653660986657467, 0.0001413130200539035, 0.01447069557673382, 0.03113945010308819, 0.0423523165735007, 0.04654843304446183, 0.04384127892295795, 0.03586297655804295, 0.02522494969228044, 0.01473594804176587, 0.006605597456080277, 0.001863641693564429, 0.000159256847770402, 0.01502795871881384, 0.0323386231293656, 0.04398329449594855, 0.04834100244236725, 0.04552959644134281, 0.0372440514963624, 0.02619635667474308, 0.01530342599496706, 0.006859977487376735, 0.001935410104444195, 0.0001653897921693557, 0.01447069557673382, 0.03113945010308819, 0.0423523165735007, 0.04654843304446183, 0.04384127892295795, 0.03586297655804295, 0.02522494969228044, 0.01473594804176587, 0.006605597456080277, 0.001863641693564429, 0.000159256847770402, 0.01284024971520857, 0.02763089812771649, 0.03758038567929434, 0.04130371625698049, 0.03890158328739785, 0.03182221421866715, 0.02238279779882984, 0.01307561558760397, 0.005861329913579828, 0.001653660986657467, 0.0001413130200539035, 0.01025761916059888, 0.02207334252415016, 0.03002163452865245, 0.03299607100161919, 0.03107709234297891, 0.02542163599166264, 0.01788082168660207, 0.01044564459125499, 0.004682408158847183, 0.001321050991849573, 0.0001128899495178914, 0.006914778815286163, 0.01487989355803276, 0.02023792843044774, 0.02224303019808677, 0.02094942465785366, 0.0171370166169049, 0.01205366713879896, 0.007041528916287183, 0.003156465085552, 0.0008905356369090305, 7.610041074477408e-05, 0.003065254786336921, 0.006596113363469354, 0.008971278567846241, 0.009860120851096311, 0.009286677986218874, 0.007596674255491598, 0.005343274438285346, 0.003121441884166165, 0.001399230542270533, 0.0003947658625615832, 3.37345784310486e-05 } };
		default:
			throw std::runtime_error("not supported order");
		}
	}
	case Figure::quadrilateral: {
		switch (integrand_order)
		{
		case 0:
		case 1:		return { { { 0, 0 } }, { 4 } };
		case 2:
		case 3:		return { { { -0.577350269189626, -0.577350269189626 }, { -0.577350269189626, 0.577350269189626 }, { 0.577350269189626, -0.577350269189626 }, { 0.577350269189626, 0.577350269189626 } }, { 1, 1, 1, 1 } };
		case 4:
		case 5:		return { { { -0.774596669241483, -0.774596669241483 }, { -0.774596669241483, 0 }, { -0.774596669241483, 0.774596669241483 }, { 0, -0.774596669241483 }, { 0, 0 }, { 0, 0.774596669241483 }, { 0.774596669241483, -0.774596669241483 }, { 0.774596669241483, 0 }, { 0.774596669241483, 0.774596669241483 } }, { 0.3086419753086403, 0.4938271604938259, 0.3086419753086403, 0.4938271604938259, 0.7901234567901235, 0.4938271604938259, 0.3086419753086403, 0.4938271604938259, 0.3086419753086403 } };
		case 6:
		case 7:		return { { { -0.861136311594052, -0.861136311594052 }, { -0.861136311594052, -0.339981043584856 }, { -0.861136311594052, 0.339981043584856 }, { -0.861136311594052, 0.861136311594052 }, { -0.339981043584856, -0.861136311594052 }, { -0.339981043584856, -0.339981043584856 }, { -0.339981043584856, 0.339981043584856 }, { -0.339981043584856, 0.861136311594052 }, { 0.339981043584856, -0.861136311594052 }, { 0.339981043584856, -0.339981043584856 }, { 0.339981043584856, 0.339981043584856 }, { 0.339981043584856, 0.861136311594052 }, { 0.861136311594052, -0.861136311594052 }, { 0.861136311594052, -0.339981043584856 }, { 0.861136311594052, 0.339981043584856 }, { 0.861136311594052, 0.861136311594052 } }, { 0.1210029932856021, 0.2268518518518519, 0.2268518518518519, 0.1210029932856021, 0.2268518518518519, 0.4252933030106941, 0.4252933030106941, 0.2268518518518519, 0.2268518518518519, 0.4252933030106941, 0.4252933030106941, 0.2268518518518519, 0.1210029932856021, 0.2268518518518519, 0.2268518518518519, 0.1210029932856021 } };
		case 8:
		case 9:		return { { { -0.906179845938664, -0.906179845938664 }, { -0.906179845938664, -0.538469310105683 }, { -0.906179845938664, 0 }, { -0.906179845938664, 0.538469310105683 }, { -0.906179845938664, 0.906179845938664 }, { -0.538469310105683, -0.906179845938664 }, { -0.538469310105683, -0.538469310105683 }, { -0.538469310105683, 0 }, { -0.538469310105683, 0.538469310105683 }, { -0.538469310105683, 0.906179845938664 }, { 0, -0.906179845938664 }, { 0, -0.538469310105683 }, { 0, 0 }, { 0, 0.538469310105683 }, { 0, 0.906179845938664 }, { 0.538469310105683, -0.906179845938664 }, { 0.538469310105683, -0.538469310105683 }, { 0.538469310105683, 0 }, { 0.538469310105683, 0.538469310105683 }, { 0.538469310105683, 0.906179845938664 }, { 0.906179845938664, -0.906179845938664 }, { 0.906179845938664, -0.538469310105683 }, { 0.906179845938664, 0 }, { 0.906179845938664, 0.538469310105683 }, { 0.906179845938664, 0.906179845938664 } }, { 0.05613434886242859, 0.1133999999999998, 0.1347850723875209, 0.1133999999999998, 0.05613434886242859, 0.1133999999999998, 0.2290854042239907, 0.2722865325507505, 0.2290854042239907, 0.1133999999999998, 0.1347850723875209, 0.2722865325507505, 0.3236345679012347, 0.2722865325507505, 0.1347850723875209, 0.1133999999999998, 0.2290854042239907, 0.2722865325507505, 0.2290854042239907, 0.1133999999999998, 0.05613434886242859, 0.1133999999999998, 0.1347850723875209, 0.1133999999999998, 0.05613434886242859 } };
		case 10:
		case 11:	return { { { -0.9324695142031521, -0.9324695142031521 }, { -0.9324695142031521, -0.661209386466264 }, { -0.9324695142031521, -0.238619186083197 }, { -0.9324695142031521, 0.238619186083197 }, { -0.9324695142031521, 0.661209386466264 }, { -0.9324695142031521, 0.9324695142031521 }, { -0.661209386466264, -0.9324695142031521 }, { -0.661209386466264, -0.661209386466264 }, { -0.661209386466264, -0.238619186083197 }, { -0.661209386466264, 0.238619186083197 }, { -0.661209386466264, 0.661209386466264 }, { -0.661209386466264, 0.9324695142031521 }, { -0.238619186083197, -0.9324695142031521 }, { -0.238619186083197, -0.661209386466264 }, { -0.238619186083197, -0.238619186083197 }, { -0.238619186083197, 0.238619186083197 }, { -0.238619186083197, 0.661209386466264 }, { -0.238619186083197, 0.9324695142031521 }, { 0.238619186083197, -0.9324695142031521 }, { 0.238619186083197, -0.661209386466264 }, { 0.238619186083197, -0.238619186083197 }, { 0.238619186083197, 0.238619186083197 }, { 0.238619186083197, 0.661209386466264 }, { 0.238619186083197, 0.9324695142031521 }, { 0.661209386466264, -0.9324695142031521 }, { 0.661209386466264, -0.661209386466264 }, { 0.661209386466264, -0.238619186083197 }, { 0.661209386466264, 0.238619186083197 }, { 0.661209386466264, 0.661209386466264 }, { 0.661209386466264, 0.9324695142031521 }, { 0.9324695142031521, -0.9324695142031521 }, { 0.9324695142031521, -0.661209386466264 }, { 0.9324695142031521, -0.238619186083197 }, { 0.9324695142031521, 0.238619186083197 }, { 0.9324695142031521, 0.661209386466264 }, { 0.9324695142031521, 0.9324695142031521 } }, { 0.02935208168898062, 0.06180729337238363, 0.08016511731780691, 0.08016511731780691, 0.06180729337238363, 0.02935208168898062, 0.06180729337238363, 0.1301489125881677, 0.168805367087588, 0.168805367087588, 0.1301489125881677, 0.06180729337238363, 0.08016511731780691, 0.168805367087588, 0.2189434501672965, 0.2189434501672965, 0.168805367087588, 0.08016511731780691, 0.08016511731780691, 0.168805367087588, 0.2189434501672965, 0.2189434501672965, 0.168805367087588, 0.08016511731780691, 0.06180729337238363, 0.1301489125881677, 0.168805367087588, 0.168805367087588, 0.1301489125881677, 0.06180729337238363, 0.02935208168898062, 0.06180729337238363, 0.08016511731780691, 0.08016511731780691, 0.06180729337238363, 0.02935208168898062 } };
		case 12:
		case 13:	return { { { -0.949107912342758, -0.949107912342758 }, { -0.949107912342758, -0.741531185599394 }, { -0.949107912342758, -0.405845151377397 }, { -0.949107912342758, 0 }, { -0.949107912342758, 0.405845151377397 }, { -0.949107912342758, 0.741531185599394 }, { -0.949107912342758, 0.949107912342758 }, { -0.741531185599394, -0.949107912342758 }, { -0.741531185599394, -0.741531185599394 }, { -0.741531185599394, -0.405845151377397 }, { -0.741531185599394, 0 }, { -0.741531185599394, 0.405845151377397 }, { -0.741531185599394, 0.741531185599394 }, { -0.741531185599394, 0.949107912342758 }, { -0.405845151377397, -0.949107912342758 }, { -0.405845151377397, -0.741531185599394 }, { -0.405845151377397, -0.405845151377397 }, { -0.405845151377397, 0 }, { -0.405845151377397, 0.405845151377397 }, { -0.405845151377397, 0.741531185599394 }, { -0.405845151377397, 0.949107912342758 }, { 0, -0.949107912342758 }, { 0, -0.741531185599394 }, { 0, -0.405845151377397 }, { 0, 0 }, { 0, 0.405845151377397 }, { 0, 0.741531185599394 }, { 0, 0.949107912342758 }, { 0.405845151377397, -0.949107912342758 }, { 0.405845151377397, -0.741531185599394 }, { 0.405845151377397, -0.405845151377397 }, { 0.405845151377397, 0 }, { 0.405845151377397, 0.405845151377397 }, { 0.405845151377397, 0.741531185599394 }, { 0.405845151377397, 0.949107912342758 }, { 0.741531185599394, -0.949107912342758 }, { 0.741531185599394, -0.741531185599394 }, { 0.741531185599394, -0.405845151377397 }, { 0.741531185599394, 0 }, { 0.741531185599394, 0.405845151377397 }, { 0.741531185599394, 0.741531185599394 }, { 0.741531185599394, 0.949107912342758 }, { 0.949107912342758, -0.949107912342758 }, { 0.949107912342758, -0.741531185599394 }, { 0.949107912342758, -0.405845151377397 }, { 0.949107912342758, 0 }, { 0.949107912342758, 0.405845151377397 }, { 0.949107912342758, 0.741531185599394 }, { 0.949107912342758, 0.949107912342758 } }, { 0.01676635646375341, 0.03621764315423957, 0.04944125117191326, 0.05411943075792766, 0.04944125117191326, 0.03621764315423957, 0.01676635646375341, 0.03621764315423957, 0.0782351060281697, 0.1067999237589047, 0.1169054370959263, 0.1067999237589047, 0.0782351060281697, 0.03621764315423957, 0.04944125117191326, 0.1067999237589047, 0.1457941874687417, 0.159589376211119, 0.1457941874687417, 0.1067999237589047, 0.04944125117191326, 0.05411943075792766, 0.1169054370959263, 0.159589376211119, 0.1746898792169926, 0.159589376211119, 0.1169054370959263, 0.05411943075792766, 0.04944125117191326, 0.1067999237589047, 0.1457941874687417, 0.159589376211119, 0.1457941874687417, 0.1067999237589047, 0.04944125117191326, 0.03621764315423957, 0.0782351060281697, 0.1067999237589047, 0.1169054370959263, 0.1067999237589047, 0.0782351060281697, 0.03621764315423957, 0.01676635646375341, 0.03621764315423957, 0.04944125117191326, 0.05411943075792766, 0.04944125117191326, 0.03621764315423957, 0.01676635646375341 } };
		case 14:
		case 15:	return { { { -0.960289856497536, -0.960289856497536 }, { -0.960289856497536, -0.796666477413627 }, { -0.960289856497536, -0.525532409916329 }, { -0.960289856497536, -0.18343464249565 }, { -0.960289856497536, 0.18343464249565 }, { -0.960289856497536, 0.525532409916329 }, { -0.960289856497536, 0.796666477413627 }, { -0.960289856497536, 0.960289856497536 }, { -0.796666477413627, -0.960289856497536 }, { -0.796666477413627, -0.796666477413627 }, { -0.796666477413627, -0.525532409916329 }, { -0.796666477413627, -0.18343464249565 }, { -0.796666477413627, 0.18343464249565 }, { -0.796666477413627, 0.525532409916329 }, { -0.796666477413627, 0.796666477413627 }, { -0.796666477413627, 0.960289856497536 }, { -0.525532409916329, -0.960289856497536 }, { -0.525532409916329, -0.796666477413627 }, { -0.525532409916329, -0.525532409916329 }, { -0.525532409916329, -0.18343464249565 }, { -0.525532409916329, 0.18343464249565 }, { -0.525532409916329, 0.525532409916329 }, { -0.525532409916329, 0.796666477413627 }, { -0.525532409916329, 0.960289856497536 }, { -0.18343464249565, -0.960289856497536 }, { -0.18343464249565, -0.796666477413627 }, { -0.18343464249565, -0.525532409916329 }, { -0.18343464249565, -0.18343464249565 }, { -0.18343464249565, 0.18343464249565 }, { -0.18343464249565, 0.525532409916329 }, { -0.18343464249565, 0.796666477413627 }, { -0.18343464249565, 0.960289856497536 }, { 0.18343464249565, -0.960289856497536 }, { 0.18343464249565, -0.796666477413627 }, { 0.18343464249565, -0.525532409916329 }, { 0.18343464249565, -0.18343464249565 }, { 0.18343464249565, 0.18343464249565 }, { 0.18343464249565, 0.525532409916329 }, { 0.18343464249565, 0.796666477413627 }, { 0.18343464249565, 0.960289856497536 }, { 0.525532409916329, -0.960289856497536 }, { 0.525532409916329, -0.796666477413627 }, { 0.525532409916329, -0.525532409916329 }, { 0.525532409916329, -0.18343464249565 }, { 0.525532409916329, 0.18343464249565 }, { 0.525532409916329, 0.525532409916329 }, { 0.525532409916329, 0.796666477413627 }, { 0.525532409916329, 0.960289856497536 }, { 0.796666477413627, -0.960289856497536 }, { 0.796666477413627, -0.796666477413627 }, { 0.796666477413627, -0.525532409916329 }, { 0.796666477413627, -0.18343464249565 }, { 0.796666477413627, 0.18343464249565 }, { 0.796666477413627, 0.525532409916329 }, { 0.796666477413627, 0.796666477413627 }, { 0.796666477413627, 0.960289856497536 }, { 0.960289856497536, -0.960289856497536 }, { 0.960289856497536, -0.796666477413627 }, { 0.960289856497536, -0.525532409916329 }, { 0.960289856497536, -0.18343464249565 }, { 0.960289856497536, 0.18343464249565 }, { 0.960289856497536, 0.525532409916329 }, { 0.960289856497536, 0.796666477413627 }, { 0.960289856497536, 0.960289856497536 } }, { 0.01024721655949217, 0.02251130661645495, 0.03175606458678213, 0.03671394852764775, 0.03671394852764775, 0.03175606458678213, 0.02251130661645495, 0.01024721655949217, 0.02251130661645495, 0.04945332448455272, 0.06976240842522279, 0.08065399492714355, 0.08065399492714355, 0.06976240842522279, 0.04945332448455272, 0.02251130661645495, 0.03175606458678213, 0.06976240842522279, 0.09841185966795399, 0.1137763131979281, 0.1137763131979281, 0.09841185966795399, 0.06976240842522279, 0.03175606458678213, 0.03671394852764775, 0.08065399492714355, 0.1137763131979281, 0.1315395267256426, 0.1315395267256426, 0.1137763131979281, 0.08065399492714355, 0.03671394852764775, 0.03671394852764775, 0.08065399492714355, 0.1137763131979281, 0.1315395267256426, 0.1315395267256426, 0.1137763131979281, 0.08065399492714355, 0.03671394852764775, 0.03175606458678213, 0.06976240842522279, 0.09841185966795399, 0.1137763131979281, 0.1137763131979281, 0.09841185966795399, 0.06976240842522279, 0.03175606458678213, 0.02251130661645495, 0.04945332448455272, 0.06976240842522279, 0.08065399492714355, 0.08065399492714355, 0.06976240842522279, 0.04945332448455272, 0.02251130661645495, 0.01024721655949217, 0.02251130661645495, 0.03175606458678213, 0.03671394852764775, 0.03671394852764775, 0.03175606458678213, 0.02251130661645495, 0.01024721655949217 } };
		case 16:
		case 17:	return { { { -0.968160239507626, -0.968160239507626 }, { -0.968160239507626, -0.836031107326636 }, { -0.968160239507626, -0.61337143270059 }, { -0.968160239507626, -0.324253423403809 }, { -0.968160239507626, 0 }, { -0.968160239507626, 0.324253423403809 }, { -0.968160239507626, 0.61337143270059 }, { -0.968160239507626, 0.836031107326636 }, { -0.968160239507626, 0.968160239507626 }, { -0.836031107326636, -0.968160239507626 }, { -0.836031107326636, -0.836031107326636 }, { -0.836031107326636, -0.61337143270059 }, { -0.836031107326636, -0.324253423403809 }, { -0.836031107326636, 0 }, { -0.836031107326636, 0.324253423403809 }, { -0.836031107326636, 0.61337143270059 }, { -0.836031107326636, 0.836031107326636 }, { -0.836031107326636, 0.968160239507626 }, { -0.61337143270059, -0.968160239507626 }, { -0.61337143270059, -0.836031107326636 }, { -0.61337143270059, -0.61337143270059 }, { -0.61337143270059, -0.324253423403809 }, { -0.61337143270059, 0 }, { -0.61337143270059, 0.324253423403809 }, { -0.61337143270059, 0.61337143270059 }, { -0.61337143270059, 0.836031107326636 }, { -0.61337143270059, 0.968160239507626 }, { -0.324253423403809, -0.968160239507626 }, { -0.324253423403809, -0.836031107326636 }, { -0.324253423403809, -0.61337143270059 }, { -0.324253423403809, -0.324253423403809 }, { -0.324253423403809, 0 }, { -0.324253423403809, 0.324253423403809 }, { -0.324253423403809, 0.61337143270059 }, { -0.324253423403809, 0.836031107326636 }, { -0.324253423403809, 0.968160239507626 }, { 0, -0.968160239507626 }, { 0, -0.836031107326636 }, { 0, -0.61337143270059 }, { 0, -0.324253423403809 }, { 0, 0 }, { 0, 0.324253423403809 }, { 0, 0.61337143270059 }, { 0, 0.836031107326636 }, { 0, 0.968160239507626 }, { 0.324253423403809, -0.968160239507626 }, { 0.324253423403809, -0.836031107326636 }, { 0.324253423403809, -0.61337143270059 }, { 0.324253423403809, -0.324253423403809 }, { 0.324253423403809, 0 }, { 0.324253423403809, 0.324253423403809 }, { 0.324253423403809, 0.61337143270059 }, { 0.324253423403809, 0.836031107326636 }, { 0.324253423403809, 0.968160239507626 }, { 0.61337143270059, -0.968160239507626 }, { 0.61337143270059, -0.836031107326636 }, { 0.61337143270059, -0.61337143270059 }, { 0.61337143270059, -0.324253423403809 }, { 0.61337143270059, 0 }, { 0.61337143270059, 0.324253423403809 }, { 0.61337143270059, 0.61337143270059 }, { 0.61337143270059, 0.836031107326636 }, { 0.61337143270059, 0.968160239507626 }, { 0.836031107326636, -0.968160239507626 }, { 0.836031107326636, -0.836031107326636 }, { 0.836031107326636, -0.61337143270059 }, { 0.836031107326636, -0.324253423403809 }, { 0.836031107326636, 0 }, { 0.836031107326636, 0.324253423403809 }, { 0.836031107326636, 0.61337143270059 }, { 0.836031107326636, 0.836031107326636 }, { 0.836031107326636, 0.968160239507626 }, { 0.968160239507626, -0.968160239507626 }, { 0.968160239507626, -0.836031107326636 }, { 0.968160239507626, -0.61337143270059 }, { 0.968160239507626, -0.324253423403809 }, { 0.968160239507626, 0 }, { 0.968160239507626, 0.324253423403809 }, { 0.968160239507626, 0.61337143270059 }, { 0.968160239507626, 0.836031107326636 }, { 0.968160239507626, 0.968160239507626 } }, { 0.006605526203548117, 0.01468206876911802, 0.02118097495063273, 0.02538581764295198, 0.02684000159064844, 0.02538581764295198, 0.02118097495063273, 0.01468206876911802, 0.006605526203548117, 0.01468206876911802, 0.03263375796243488, 0.04707884296259617, 0.05642492496569134, 0.05965713207003355, 0.05642492496569134, 0.04707884296259617, 0.03263375796243488, 0.01468206876911802, 0.02118097495063273, 0.04707884296259617, 0.06791793507962328, 0.08140098926681667, 0.08606390828653478, 0.08140098926681667, 0.06791793507962328, 0.04707884296259617, 0.02118097495063273, 0.02538581764295198, 0.05642492496569134, 0.08140098926681667, 0.09756069653543355, 0.1031492972582194, 0.09756069653543355, 0.08140098926681667, 0.05642492496569134, 0.02538581764295198, 0.02684000159064844, 0.05965713207003355, 0.08606390828653478, 0.1031492972582194, 0.1090580315916482, 0.1031492972582194, 0.08606390828653478, 0.05965713207003355, 0.02684000159064844, 0.02538581764295198, 0.05642492496569134, 0.08140098926681667, 0.09756069653543355, 0.1031492972582194, 0.09756069653543355, 0.08140098926681667, 0.05642492496569134, 0.02538581764295198, 0.02118097495063273, 0.04707884296259617, 0.06791793507962328, 0.08140098926681667, 0.08606390828653478, 0.08140098926681667, 0.06791793507962328, 0.04707884296259617, 0.02118097495063273, 0.01468206876911802, 0.03263375796243488, 0.04707884296259617, 0.05642492496569134, 0.05965713207003355, 0.05642492496569134, 0.04707884296259617, 0.03263375796243488, 0.01468206876911802, 0.006605526203548117, 0.01468206876911802, 0.02118097495063273, 0.02538581764295198, 0.02684000159064844, 0.02538581764295198, 0.02118097495063273, 0.01468206876911802, 0.006605526203548117 } };
		case 18:
		case 19:	return { { { -0.973906528517172, -0.973906528517172 }, { -0.973906528517172, -0.865063366688985 }, { -0.973906528517172, -0.679409568299024 }, { -0.973906528517172, -0.433395394129247 }, { -0.973906528517172, -0.148874338981631 }, { -0.973906528517172, 0.148874338981631 }, { -0.973906528517172, 0.433395394129247 }, { -0.973906528517172, 0.679409568299024 }, { -0.973906528517172, 0.865063366688985 }, { -0.973906528517172, 0.973906528517172 }, { -0.865063366688985, -0.973906528517172 }, { -0.865063366688985, -0.865063366688985 }, { -0.865063366688985, -0.679409568299024 }, { -0.865063366688985, -0.433395394129247 }, { -0.865063366688985, -0.148874338981631 }, { -0.865063366688985, 0.148874338981631 }, { -0.865063366688985, 0.433395394129247 }, { -0.865063366688985, 0.679409568299024 }, { -0.865063366688985, 0.865063366688985 }, { -0.865063366688985, 0.973906528517172 }, { -0.679409568299024, -0.973906528517172 }, { -0.679409568299024, -0.865063366688985 }, { -0.679409568299024, -0.679409568299024 }, { -0.679409568299024, -0.433395394129247 }, { -0.679409568299024, -0.148874338981631 }, { -0.679409568299024, 0.148874338981631 }, { -0.679409568299024, 0.433395394129247 }, { -0.679409568299024, 0.679409568299024 }, { -0.679409568299024, 0.865063366688985 }, { -0.679409568299024, 0.973906528517172 }, { -0.433395394129247, -0.973906528517172 }, { -0.433395394129247, -0.865063366688985 }, { -0.433395394129247, -0.679409568299024 }, { -0.433395394129247, -0.433395394129247 }, { -0.433395394129247, -0.148874338981631 }, { -0.433395394129247, 0.148874338981631 }, { -0.433395394129247, 0.433395394129247 }, { -0.433395394129247, 0.679409568299024 }, { -0.433395394129247, 0.865063366688985 }, { -0.433395394129247, 0.973906528517172 }, { -0.148874338981631, -0.973906528517172 }, { -0.148874338981631, -0.865063366688985 }, { -0.148874338981631, -0.679409568299024 }, { -0.148874338981631, -0.433395394129247 }, { -0.148874338981631, -0.148874338981631 }, { -0.148874338981631, 0.148874338981631 }, { -0.148874338981631, 0.433395394129247 }, { -0.148874338981631, 0.679409568299024 }, { -0.148874338981631, 0.865063366688985 }, { -0.148874338981631, 0.973906528517172 }, { 0.148874338981631, -0.973906528517172 }, { 0.148874338981631, -0.865063366688985 }, { 0.148874338981631, -0.679409568299024 }, { 0.148874338981631, -0.433395394129247 }, { 0.148874338981631, -0.148874338981631 }, { 0.148874338981631, 0.148874338981631 }, { 0.148874338981631, 0.433395394129247 }, { 0.148874338981631, 0.679409568299024 }, { 0.148874338981631, 0.865063366688985 }, { 0.148874338981631, 0.973906528517172 }, { 0.433395394129247, -0.973906528517172 }, { 0.433395394129247, -0.865063366688985 }, { 0.433395394129247, -0.679409568299024 }, { 0.433395394129247, -0.433395394129247 }, { 0.433395394129247, -0.148874338981631 }, { 0.433395394129247, 0.148874338981631 }, { 0.433395394129247, 0.433395394129247 }, { 0.433395394129247, 0.679409568299024 }, { 0.433395394129247, 0.865063366688985 }, { 0.433395394129247, 0.973906528517172 }, { 0.679409568299024, -0.973906528517172 }, { 0.679409568299024, -0.865063366688985 }, { 0.679409568299024, -0.679409568299024 }, { 0.679409568299024, -0.433395394129247 }, { 0.679409568299024, -0.148874338981631 }, { 0.679409568299024, 0.148874338981631 }, { 0.679409568299024, 0.433395394129247 }, { 0.679409568299024, 0.679409568299024 }, { 0.679409568299024, 0.865063366688985 }, { 0.679409568299024, 0.973906528517172 }, { 0.865063366688985, -0.973906528517172 }, { 0.865063366688985, -0.865063366688985 }, { 0.865063366688985, -0.679409568299024 }, { 0.865063366688985, -0.433395394129247 }, { 0.865063366688985, -0.148874338981631 }, { 0.865063366688985, 0.148874338981631 }, { 0.865063366688985, 0.433395394129247 }, { 0.865063366688985, 0.679409568299024 }, { 0.865063366688985, 0.865063366688985 }, { 0.865063366688985, 0.973906528517172 }, { 0.973906528517172, -0.973906528517172 }, { 0.973906528517172, -0.865063366688985 }, { 0.973906528517172, -0.679409568299024 }, { 0.973906528517172, -0.433395394129247 }, { 0.973906528517172, -0.148874338981631 }, { 0.973906528517172, 0.148874338981631 }, { 0.973906528517172, 0.433395394129247 }, { 0.973906528517172, 0.679409568299024 }, { 0.973906528517172, 0.865063366688985 }, { 0.973906528517172, 0.973906528517172 } }, { 0.004445068151927624, 0.009964122356616333, 0.01460678230864107, 0.01795237415398759, 0.01970299733751538, 0.01970299733751538, 0.01795237415398759, 0.01460678230864107, 0.009964122356616333, 0.004445068151927624, 0.009964122356616333, 0.02233570576292887, 0.03274275245850679, 0.04024227448222971, 0.04416649409029931, 0.04416649409029931, 0.04024227448222971, 0.03274275245850679, 0.02233570576292887, 0.009964122356616333, 0.01460678230864107, 0.03274275245850679, 0.04799883424048429, 0.05899266608023896, 0.0647453274281109, 0.0647453274281109, 0.05899266608023896, 0.04799883424048429, 0.03274275245850679, 0.01460678230864107, 0.01795237415398759, 0.04024227448222971, 0.05899266608023896, 0.07250456612796818, 0.07957483846557158, 0.07957483846557158, 0.07250456612796818, 0.05899266608023896, 0.04024227448222971, 0.01795237415398759, 0.01970299733751538, 0.04416649409029931, 0.0647453274281109, 0.07957483846557158, 0.08733456739325582, 0.08733456739325582, 0.07957483846557158, 0.0647453274281109, 0.04416649409029931, 0.01970299733751538, 0.01970299733751538, 0.04416649409029931, 0.0647453274281109, 0.07957483846557158, 0.08733456739325582, 0.08733456739325582, 0.07957483846557158, 0.0647453274281109, 0.04416649409029931, 0.01970299733751538, 0.01795237415398759, 0.04024227448222971, 0.05899266608023896, 0.07250456612796818, 0.07957483846557158, 0.07957483846557158, 0.07250456612796818, 0.05899266608023896, 0.04024227448222971, 0.01795237415398759, 0.01460678230864107, 0.03274275245850679, 0.04799883424048429, 0.05899266608023896, 0.0647453274281109, 0.0647453274281109, 0.05899266608023896, 0.04799883424048429, 0.03274275245850679, 0.01460678230864107, 0.009964122356616333, 0.02233570576292887, 0.03274275245850679, 0.04024227448222971, 0.04416649409029931, 0.04416649409029931, 0.04024227448222971, 0.03274275245850679, 0.02233570576292887, 0.009964122356616333, 0.004445068151927624, 0.009964122356616333, 0.01460678230864107, 0.01795237415398759, 0.01970299733751538, 0.01970299733751538, 0.01795237415398759, 0.01460678230864107, 0.009964122356616333, 0.004445068151927624 } };
		case 20:
		case 21:	return { { { -0.978228658146057, -0.978228658146057 }, { -0.978228658146057, -0.887062599768095 }, { -0.978228658146057, -0.730152005574049 }, { -0.978228658146057, -0.519096129206812 }, { -0.978228658146057, -0.269543155952345 }, { -0.978228658146057, 0 }, { -0.978228658146057, 0.269543155952345 }, { -0.978228658146057, 0.519096129206812 }, { -0.978228658146057, 0.730152005574049 }, { -0.978228658146057, 0.887062599768095 }, { -0.978228658146057, 0.978228658146057 }, { -0.887062599768095, -0.978228658146057 }, { -0.887062599768095, -0.887062599768095 }, { -0.887062599768095, -0.730152005574049 }, { -0.887062599768095, -0.519096129206812 }, { -0.887062599768095, -0.269543155952345 }, { -0.887062599768095, 0 }, { -0.887062599768095, 0.269543155952345 }, { -0.887062599768095, 0.519096129206812 }, { -0.887062599768095, 0.730152005574049 }, { -0.887062599768095, 0.887062599768095 }, { -0.887062599768095, 0.978228658146057 }, { -0.730152005574049, -0.978228658146057 }, { -0.730152005574049, -0.887062599768095 }, { -0.730152005574049, -0.730152005574049 }, { -0.730152005574049, -0.519096129206812 }, { -0.730152005574049, -0.269543155952345 }, { -0.730152005574049, 0 }, { -0.730152005574049, 0.269543155952345 }, { -0.730152005574049, 0.519096129206812 }, { -0.730152005574049, 0.730152005574049 }, { -0.730152005574049, 0.887062599768095 }, { -0.730152005574049, 0.978228658146057 }, { -0.519096129206812, -0.978228658146057 }, { -0.519096129206812, -0.887062599768095 }, { -0.519096129206812, -0.730152005574049 }, { -0.519096129206812, -0.519096129206812 }, { -0.519096129206812, -0.269543155952345 }, { -0.519096129206812, 0 }, { -0.519096129206812, 0.269543155952345 }, { -0.519096129206812, 0.519096129206812 }, { -0.519096129206812, 0.730152005574049 }, { -0.519096129206812, 0.887062599768095 }, { -0.519096129206812, 0.978228658146057 }, { -0.269543155952345, -0.978228658146057 }, { -0.269543155952345, -0.887062599768095 }, { -0.269543155952345, -0.730152005574049 }, { -0.269543155952345, -0.519096129206812 }, { -0.269543155952345, -0.269543155952345 }, { -0.269543155952345, 0 }, { -0.269543155952345, 0.269543155952345 }, { -0.269543155952345, 0.519096129206812 }, { -0.269543155952345, 0.730152005574049 }, { -0.269543155952345, 0.887062599768095 }, { -0.269543155952345, 0.978228658146057 }, { 0, -0.978228658146057 }, { 0, -0.887062599768095 }, { 0, -0.730152005574049 }, { 0, -0.519096129206812 }, { 0, -0.269543155952345 }, { 0, 0 }, { 0, 0.269543155952345 }, { 0, 0.519096129206812 }, { 0, 0.730152005574049 }, { 0, 0.887062599768095 }, { 0, 0.978228658146057 }, { 0.269543155952345, -0.978228658146057 }, { 0.269543155952345, -0.887062599768095 }, { 0.269543155952345, -0.730152005574049 }, { 0.269543155952345, -0.519096129206812 }, { 0.269543155952345, -0.269543155952345 }, { 0.269543155952345, 0 }, { 0.269543155952345, 0.269543155952345 }, { 0.269543155952345, 0.519096129206812 }, { 0.269543155952345, 0.730152005574049 }, { 0.269543155952345, 0.887062599768095 }, { 0.269543155952345, 0.978228658146057 }, { 0.519096129206812, -0.978228658146057 }, { 0.519096129206812, -0.887062599768095 }, { 0.519096129206812, -0.730152005574049 }, { 0.519096129206812, -0.519096129206812 }, { 0.519096129206812, -0.269543155952345 }, { 0.519096129206812, 0 }, { 0.519096129206812, 0.269543155952345 }, { 0.519096129206812, 0.519096129206812 }, { 0.519096129206812, 0.730152005574049 }, { 0.519096129206812, 0.887062599768095 }, { 0.519096129206812, 0.978228658146057 }, { 0.730152005574049, -0.978228658146057 }, { 0.730152005574049, -0.887062599768095 }, { 0.730152005574049, -0.730152005574049 }, { 0.730152005574049, -0.519096129206812 }, { 0.730152005574049, -0.269543155952345 }, { 0.730152005574049, 0 }, { 0.730152005574049, 0.269543155952345 }, { 0.730152005574049, 0.519096129206812 }, { 0.730152005574049, 0.730152005574049 }, { 0.730152005574049, 0.887062599768095 }, { 0.730152005574049, 0.978228658146057 }, { 0.887062599768095, -0.978228658146057 }, { 0.887062599768095, -0.887062599768095 }, { 0.887062599768095, -0.730152005574049 }, { 0.887062599768095, -0.519096129206812 }, { 0.887062599768095, -0.269543155952345 }, { 0.887062599768095, 0 }, { 0.887062599768095, 0.269543155952345 }, { 0.887062599768095, 0.519096129206812 }, { 0.887062599768095, 0.730152005574049 }, { 0.887062599768095, 0.887062599768095 }, { 0.887062599768095, 0.978228658146057 }, { 0.978228658146057, -0.978228658146057 }, { 0.978228658146057, -0.887062599768095 }, { 0.978228658146057, -0.730152005574049 }, { 0.978228658146057, -0.519096129206812 }, { 0.978228658146057, -0.269543155952345 }, { 0.978228658146057, 0 }, { 0.978228658146057, 0.269543155952345 }, { 0.978228658146057, 0.519096129206812 }, { 0.978228658146057, 0.730152005574049 }, { 0.978228658146057, 0.887062599768095 }, { 0.978228658146057, 0.978228658146057 } }, { 0.00309898936476797, 0.006990879226030937, 0.01037050911011677, 0.01298156273526248, 0.01462995242450422, 0.0151933485109832, 0.01462995242450422, 0.01298156273526248, 0.01037050911011677, 0.006990879226030937, 0.00309898936476797, 0.006990879226030937, 0.01577042919494179, 0.02339439351599974, 0.02928455911437395, 0.03300309179665262, 0.03427403323380979, 0.03300309179665262, 0.02928455911437395, 0.02339439351599974, 0.01577042919494179, 0.006990879226030937, 0.01037050911011677, 0.02339439351599974, 0.03470404268749963, 0.04344171559287417, 0.04895791402958097, 0.05084327198332528, 0.04895791402958097, 0.04344171559287417, 0.03470404268749963, 0.02339439351599974, 0.01037050911011677, 0.01298156273526248, 0.02928455911437395, 0.04344171559287417, 0.05437933184458445, 0.0612843810862277, 0.0636444284373343, 0.0612843810862277, 0.05437933184458445, 0.04344171559287417, 0.02928455911437395, 0.01298156273526248, 0.01462995242450422, 0.03300309179665262, 0.04895791402958097, 0.0612843810862277, 0.0690662286152384, 0.0717259531160859, 0.0690662286152384, 0.0612843810862277, 0.04895791402958097, 0.03300309179665262, 0.01462995242450422, 0.0151933485109832, 0.03427403323380979, 0.05084327198332528, 0.0636444284373343, 0.0717259531160859, 0.07448810299272479, 0.0717259531160859, 0.0636444284373343, 0.05084327198332528, 0.03427403323380979, 0.0151933485109832, 0.01462995242450422, 0.03300309179665262, 0.04895791402958097, 0.0612843810862277, 0.0690662286152384, 0.0717259531160859, 0.0690662286152384, 0.0612843810862277, 0.04895791402958097, 0.03300309179665262, 0.01462995242450422, 0.01298156273526248, 0.02928455911437395, 0.04344171559287417, 0.05437933184458445, 0.0612843810862277, 0.0636444284373343, 0.0612843810862277, 0.05437933184458445, 0.04344171559287417, 0.02928455911437395, 0.01298156273526248, 0.01037050911011677, 0.02339439351599974, 0.03470404268749963, 0.04344171559287417, 0.04895791402958097, 0.05084327198332528, 0.04895791402958097, 0.04344171559287417, 0.03470404268749963, 0.02339439351599974, 0.01037050911011677, 0.006990879226030937, 0.01577042919494179, 0.02339439351599974, 0.02928455911437395, 0.03300309179665262, 0.03427403323380979, 0.03300309179665262, 0.02928455911437395, 0.02339439351599974, 0.01577042919494179, 0.006990879226030937, 0.00309898936476797, 0.006990879226030937, 0.01037050911011677, 0.01298156273526248, 0.01462995242450422, 0.0151933485109832, 0.01462995242450422, 0.01298156273526248, 0.01037050911011677, 0.006990879226030937, 0.00309898936476797 } };
		default:
			throw std::runtime_error("not supported order");
		}
	}
	default:
		throw std::runtime_error("not supported figure");
		break;
	}
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<space_dimension>> ReferenceGeometry<space_dimension>::reference_post_nodes(const ushort post_order) const {
	std::vector<Space_Vector_> reference_post_nodes;

	switch (this->figure_) {
	case Figure::triangle: {
		const ushort num_reference_post_point = static_cast<ushort>((post_order + 2) * (post_order + 3) * 0.5);
		reference_post_nodes.reserve(num_reference_post_point);

		const double delta = 2.0 / (post_order + 1);

		const double X0_start_coord = -1.0;
		const double X1_start_coord = -1.0;

		for (ushort j = 0; j <= post_order + 1; ++j) {
			for (ushort i = 0; i <= post_order + 1 - j; ++i) {
				const double X0_coord = X0_start_coord + delta * i;
				const double X1_coord = X1_start_coord + delta * j;

				reference_post_nodes.push_back({ X0_coord, X1_coord });
			}
		}
		break;
	}
	case Figure::quadrilateral: {
		const ushort num_reference_post_point = (post_order + 2) * (post_order + 2);
		reference_post_nodes.reserve(num_reference_post_point);

		const double delta = 2.0 / (post_order + 1);

		const double X0_start_coord = -1.0;
		const double X1_start_coord = -1.0;

		for (ushort j = 0; j <= post_order + 1; ++j) {
			for (ushort i = 0; i <= post_order + 1; ++i) {
				const double X0_coord = X0_start_coord + delta * i;
				const double X1_coord = X1_start_coord + delta * j;

				reference_post_nodes.push_back({ X0_coord, X1_coord });
			}
		}
		break;
	}
	default:
		throw std::runtime_error("not supproted figure");
		break;
	}

	return reference_post_nodes;
}

template <ushort space_dimension>
std::vector<std::vector<ushort>> ReferenceGeometry<space_dimension>::reference_connectivity(const ushort post_order) const {
	std::vector<std::vector<ushort>> sub_element_indexes;

	switch (this->figure_) {
	case Figure::triangle: {
		const ushort num_sub_element = (post_order + 1) * (post_order + 1);
		sub_element_indexes.resize(num_sub_element);

		constexpr ushort num_sub_element_node = 3;

		ushort count = 0;
		for (ushort j = 0; j < post_order + 1; j++) {
			for (ushort i = 0; i < post_order + 1 - j; i++) {
				//	next_index		a_(j+1) + i  ------------	a_(j+1) + 1			// a_0		= 0
				//						|							|				// a_(j+1)	= a_j + (order + 2) - j
				//						|							|				// By recurrence relation
				//	start_index		a_j + i		-------------	a_j + i + 1			// a_j		= (order + 3) * j - j*(j+1)/2

				const ushort a_j = (post_order + 3) * j - static_cast<ushort>(j * (j + 1) * 0.5);
				const ushort a_jp1 = (post_order + 3) * (j + 1) - static_cast<ushort>((j + 1) * (j + 2) * 0.5);

				const ushort start_index = a_j + i;
				const ushort next_index = a_jp1 + i;

				sub_element_indexes[count].resize(num_sub_element_node);
				sub_element_indexes[count][0] = start_index;
				sub_element_indexes[count][1] = start_index + 1;
				sub_element_indexes[count][2] = next_index;
				count++;

				if (i == post_order - j)
					continue;

				sub_element_indexes[count].resize(num_sub_element_node);
				sub_element_indexes[count][0] = start_index + 1;
				sub_element_indexes[count][1] = next_index + 1;
				sub_element_indexes[count][2] = next_index;
				count++;
			}
		}

		return sub_element_indexes;
	}

	case Figure::quadrilateral: {
		const ushort num_simplex = 2 * (post_order + 1) * (post_order + 1);
		sub_element_indexes.resize(num_simplex);

		const ushort num_sub_element_node = 3;

		ushort count = 0;
		for (ushort j = 0; j <= post_order; j++) {
			for (ushort i = 0; i <= post_order; i++) {
				//	next_index		a_(j+1) + i  ------------	a_(j+1) + 1			// a_0		= 0
				//						|							|				// a_(j+1)	= a_j + order + 2
				//						|							|				// By recurrence relation
				//	start_index		a_j + i		-------------	a_j + i + 1			// a_j		= (order + 2) * j

				const ushort a_j = (post_order + 2) * j;
				const ushort a_jp1 = (post_order + 2) * (j + 1);

				const ushort start_index = a_j + i;
				const ushort next_index = a_jp1 + i;

				sub_element_indexes[count].resize(num_sub_element_node);
				sub_element_indexes[count][0] = start_index;
				sub_element_indexes[count][1] = start_index + 1;
				sub_element_indexes[count][2] = next_index;
				count++;

				sub_element_indexes[count].resize(num_sub_element_node);
				sub_element_indexes[count][0] = start_index + 1;
				sub_element_indexes[count][1] = next_index + 1;
				sub_element_indexes[count][2] = next_index;
				count++;
			}
		}

		return sub_element_indexes;
	}

	default:
		throw std::runtime_error("not supported figure");
		return sub_element_indexes;
	}
}

//template <ushort space_dimension>
//Dynamic_Vector_Function_<Polynomial<space_dimension>> operator*(const Dynamic_Matrix_& A, const Dynamic_Vector_Function_<Polynomial<space_dimension>>& v) {
//	const auto [num_row, num_column] = A.size();
//	const auto range_dimension = v.range_dimension();
//	dynamic_require(num_column == range_dimension, "number of column should be same with range dimension");
//
//	Dynamic_Vector_Function_<Polynomial<space_dimension>> result(num_row);
//	for (size_t i = 0; i < num_row; ++i)
//		for (size_t j = 0; j < num_column; ++j)
//			result[i] += A.at(i, j) * v.at(j);
//
//	return result;
//}




template <ushort space_dimension>
std::array<double, space_dimension> Geometry<space_dimension>::coordinate_projected_volume(void) const {
	if constexpr (space_dimension == 2) {
		double x_projected_volume = 0.0;
		double y_projected_volume = 0.0;

		const auto faces_nodes = this->calculate_faces_nodes();
		for (const auto& face_nodes : faces_nodes) {
			const auto& start_node = face_nodes[0];
			const auto& end_node = face_nodes[1];
			const auto node_to_node = end_node - start_node;

			x_projected_volume += std::abs(node_to_node.at(0));
			y_projected_volume += std::abs(node_to_node.at(1));
		}

		return { 0.5 * x_projected_volume, 0.5 * y_projected_volume };
	}
	else {
		throw std::runtime_error("not supproted dimension");
		return {};
	}
}

template <ushort space_dimension>
std::vector<Geometry<space_dimension>> Geometry<space_dimension>::faces_geometry(void) const {
	const auto face_reference_geometries = this->reference_geometry_.face_reference_geometries();
	const auto faces_node_index_orders = this->reference_geometry_.face_node_index_orders_set();
	const auto num_face = faces_node_index_orders.size();

	std::vector<Geometry> faces_geometry;
	faces_geometry.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i) {
		const auto& reference_geometry = face_reference_geometries[i];

		const auto& node_index_orders = faces_node_index_orders[i];
		const auto num_node = node_index_orders.size();

		std::vector<Space_Vector_> face_nodes(num_node);
		for (size_t j = 0; j < num_node; ++j) 
			face_nodes[j] = this->nodes_[node_index_orders[j]];

		faces_geometry.push_back({ reference_geometry,std::move(face_nodes) });
	}

	return faces_geometry;
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<space_dimension>> Geometry<space_dimension>::post_nodes(const ushort post_order) const {
	return this->reference_geometry_.post_nodes(this->mapping_function_, post_order);
}

template <ushort space_dimension>
std::vector<Euclidean_Vector<space_dimension>> Geometry<space_dimension>::vertex_nodes(void) const {
	const auto vertex_node_index_orders = this->reference_geometry_.vertex_node_index_orders();
	const auto num_vertex_node = vertex_node_index_orders.size();

	std::vector<Space_Vector_> vertex_nodes(num_vertex_node);
	for (size_t i = 0; i < num_vertex_node; ++i)
		vertex_nodes[i] = this->nodes_[vertex_node_index_orders[i]];

	return vertex_nodes;
}

template <ushort space_dimension>
std::vector<std::vector<typename Geometry<space_dimension>::Space_Vector_>> Geometry<space_dimension>::calculate_faces_nodes(void) const {
	const auto faces_node_index_orders = this->reference_geometry_.face_node_index_orders_set();
	const auto num_face = faces_node_index_orders.size();

	std::vector<std::vector<Space_Vector_>> faces_nodes(num_face);
	for (size_t i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = faces_node_index_orders[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_nodes = faces_nodes[i];
		face_nodes.resize(num_node);

		for (size_t j = 0; j < face_node_index_orders.size(); ++j)
			face_nodes[j] = this->nodes_[face_node_index_orders[j]];
	}

	return faces_nodes;
}

template <ushort space_dimension>
bool Geometry<space_dimension>::is_axis_parallel(const Geometry& other, const ushort axis_tag) const {
	if (this->reference_geometry_ != other.reference_geometry_)
		return false;

	if (this->nodes_.size() != other.nodes_.size())
		return false;

	for (const auto& node : other.nodes_) {
		if (this->is_axis_parallel_node(node, axis_tag))
			continue;
		else
			return false;
	}
	return true;
}

template <ushort space_dimension>
template <ushort order>
auto Geometry<space_dimension>::initial_basis_function(void) const {
	constexpr auto num_basis = ms::combination_with_repetition(1 + space_dimension, order);

	std::array<Polynomial<space_dimension>, num_basis> initial_basis_set = { 0 };
	
	if (space_dimension == 2) {
		ushort index = 0;

		//1 (x - x_c) (y - y_c)  ...
		Polynomial<space_dimension> x("x0");
		Polynomial<space_dimension> y("x1");

		const auto center_node = this->center_node();
		const auto x_c = center_node.at(0);
		const auto y_c = center_node.at(1);

		for (ushort a = 0; a <= order; ++a)
			for (ushort b = 0; b <= a; ++b)
				initial_basis_set[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
	}
	else
		throw std::runtime_error("not supported space dimension");

	Vector_Function<Polynomial<space_dimension>, num_basis> initial_basis_function = initial_basis_set;
	return initial_basis_function;
}

template <ushort space_dimension>
template <ushort order>
auto Geometry<space_dimension>::orthonormal_basis_function(void) const {
	const auto initial_basis_set = this->initial_basis_function<order>();
	return ms::Gram_Schmidt_process(initial_basis_set, *this);
}

template <ushort space_dimension>
const Quadrature_Rule<space_dimension>& Geometry<space_dimension>::get_quadrature_rule(const ushort integrand_order) const {
	if (this->integrand_order_to_quadrature_rule_.find(integrand_order) == this->integrand_order_to_quadrature_rule_.end())
		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->reference_geometry_.quadrature_rule(this->mapping_function_, integrand_order));

	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
}

template <ushort space_dimension>
bool Geometry<space_dimension>::is_axis_parallel_node(const Space_Vector_& node, const size_t axis_tag) const {
	for (const auto& my_node : this->nodes_) {
		if (my_node.is_axis_translation(node,axis_tag))
			return true;
	}
	return false;
}

template <ushort space_dimension>
Euclidean_Vector<space_dimension> Geometry<space_dimension>::center_node(void) const {
	const auto center_temp = this->mapping_function_(this->reference_geometry_.center_node());

	std::array<double, space_dimension> center;
	for (size_t i = 0; i < space_dimension; ++i)
		center[i] = center_temp.at(i);

	return center;
}

template <ushort space_dimension>
Euclidean_Vector<space_dimension> Geometry<space_dimension>::normalized_normal_vector(const Space_Vector_& node) const {
	const auto normal_vector_function = this->reference_geometry_.normal_vector_function(this->mapping_function_);
	auto normal_vector = normal_vector_function(node);

	return normal_vector.be_normalize();
}

template <ushort space_dimension>
double Geometry<space_dimension>::volume(void) const {
	const auto& quadrature_rule = this->get_quadrature_rule(0);

	auto volume = 0.0;
	for (const auto weight : quadrature_rule.weights)
		volume += weight;

	return volume;
}

template <ushort space_dimension>
ElementType Element<space_dimension>::type(void) const {
	return this->element_type_;
}

template <ushort space_dimension>
Euclidean_Vector<space_dimension> Element<space_dimension>::normalized_normal_vector(const Element& owner_cell_element, const Euclidean_Vector<space_dimension>& node) const {
	auto normal_vector = this->geometry_.normalized_normal_vector(node);
	
	const auto face_type = this->check_face_type(owner_cell_element);
	if (face_type == FaceType::inward_face)
		normal_vector *= -1;

	return normal_vector;
}

template <ushort space_dimension>
std::vector<uint> Element<space_dimension>::vertex_node_indexes(void) const {
	const auto num_vertex = this->geometry_.reference_geometry_.num_vertex();

	return { this->node_indexes_.begin(), this->node_indexes_.begin() + num_vertex };
}

template <ushort space_dimension>
std::vector<Element<space_dimension>> Element<space_dimension>::make_inner_face_elements(void) const {
	dynamic_require(this->element_type_ == ElementType::cell, "make inner face elements should be called from cell element");

	auto faces_geometry = this->geometry_.faces_geometry();
	auto faces_node_indexes = this->face_node_indexes_set();

	const auto num_face = faces_geometry.size();
	std::vector<Element<space_dimension>> inner_face_elements;
	inner_face_elements.reserve(num_face);

	for (ushort i = 0; i < num_face; ++i) 
		inner_face_elements.push_back({ ElementType::inner_face, std::move(faces_geometry[i]), std::move(faces_node_indexes[i]) });

	return inner_face_elements;
}


template <ushort space_dimension>
std::vector<std::vector<uint>> Element<space_dimension>::face_node_indexes_set(void) const {
	const auto face_node_index_orders_set = this->geometry_.reference_geometry_.face_node_index_orders_set();
	const auto num_face = face_node_index_orders_set.size();

	std::vector<std::vector<uint>> face_node_indexes_set(num_face);
	for (ushort i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = face_node_index_orders_set[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_node_indexes = face_node_indexes_set[i];
		face_node_indexes.resize(num_node);

		for (ushort j = 0; j < num_node; ++j)
			face_node_indexes[j] = this->node_indexes_[face_node_index_orders[j]];
	}

	return face_node_indexes_set;
}


template <ushort space_dimension>
std::vector<std::vector<uint>> Element<space_dimension>::face_vertex_node_indexes_set(void) const {
	const auto face_vnode_index_orders_set = this->geometry_.reference_geometry_.face_vertex_node_index_orders_set();
	const auto num_face = face_vnode_index_orders_set.size();

	std::vector<std::vector<uint>> face_vnode_indexes_set(num_face);
	for (uint i = 0; i < num_face; ++i) {
		const auto& face_node_index_orders = face_vnode_index_orders_set[i];
		const auto num_node = face_node_index_orders.size();

		auto& face_node_indexes = face_vnode_indexes_set[i];
		face_node_indexes.resize(num_node);

		for (ushort j = 0; j < num_node; ++j)
			face_node_indexes[j] = this->node_indexes_[face_node_index_orders[j]];
	}

	return face_vnode_indexes_set;
}

template <ushort space_dimension>
bool Element<space_dimension>::is_periodic_pair(const Element& other) const {
	dynamic_require(this->is_periodic_boundary() && other.is_periodic_boundary(), "both elemets should be periodic boundary");

	if (this->element_type_ != other.element_type_)
		return false;

	ushort axis_tag;
	if (this->element_type_ == ElementType::x_periodic)
		axis_tag = 0;
	else
		axis_tag = 1;

	if (this->geometry_.is_axis_parallel(other.geometry_, axis_tag))
		return true;
	else
		return false;
}

template <ushort space_dimension>
std::vector<std::pair<uint, uint>> Element<space_dimension>::find_periodic_vnode_index_pairs(const Element& other) const {
	//both element should be periodic pair
	const auto this_vnode_indexes = this->vertex_node_indexes();
	const auto other_vnode_indexes = other.vertex_node_indexes();

	const auto this_vnodes = this->geometry_.vertex_nodes();
	const auto other_vnodes = other.geometry_.vertex_nodes();

	const auto this_num_vnode = this_vnode_indexes.size();
	const auto other_num_vnode = other_vnode_indexes.size();
	dynamic_require(this_num_vnode == other_num_vnode, "periodic pair should have same number of vertex node");
	
	ushort axis_tag = 0;
	if (this->element_type_ == ElementType::x_periodic)
		axis_tag = 0;
	else if (this->element_type_ == ElementType::y_periodic)
		axis_tag = 1;
	else
		throw std::runtime_error("wrong element type");
	dynamic_require(this->element_type_ == other.element_type_, "periodic pair should have same element type");

	std::unordered_set<uint> matched_other_vnode_index;
	matched_other_vnode_index.reserve(other_num_vnode);

	std::vector<std::pair<uint, uint>> periodic_vnode_index_pairs;
	periodic_vnode_index_pairs.reserve(this_num_vnode);

	for (ushort i = 0; i < this_num_vnode; ++i) {
		const auto& this_vnode = this_vnodes[i];
		const auto this_vnode_index = this_vnode_indexes[i];
		for (ushort j = 0; j < other_num_vnode; ++j) {
			const auto& other_vnode = other_vnodes[j];
			const auto other_vnode_index = other_vnode_indexes[j];

			if (matched_other_vnode_index.find(other_vnode_index) != matched_other_vnode_index.end())
				continue;

			if (this_vnode.is_axis_translation(other_vnode, axis_tag)) {
				periodic_vnode_index_pairs.push_back(std::make_pair(this_vnode_index, other_vnode_index));
				matched_other_vnode_index.emplace(other_vnode_index);
			}
		}
	}
	
	dynamic_require(periodic_vnode_index_pairs.size() == this_num_vnode, "every vnode should have pair");
	return periodic_vnode_index_pairs;
}


template <ushort space_dimension>
bool Element<space_dimension>::is_periodic_boundary(void) const {
	switch (this->element_type_)	{
		case ElementType::x_periodic:
		case ElementType::y_periodic:	return true;
		default:						return false;
	}
}

template <ushort space_dimension>
FaceType Element<space_dimension>::check_face_type(const Element& owner_cell_element) const {
	dynamic_require(this->element_type_ != ElementType::cell, "face or boundary element should be use this method");

	const auto this_vnode_indexes = this->vertex_node_indexes();
	const auto set_of_face_vnode_indexes = owner_cell_element.face_vertex_node_indexes_set();

	for (const auto& face_vnode_indexes : set_of_face_vnode_indexes) {
		if (!std::is_permutation(face_vnode_indexes.begin(), face_vnode_indexes.end(), this_vnode_indexes.begin(), this_vnode_indexes.end()))
			continue;
		else if (this_vnode_indexes == face_vnode_indexes)
			return FaceType::inward_face;
		else
			return FaceType::outward_face;
	}

	throw std::runtime_error("this face is not input cell element face!");
	return FaceType::not_my_face;
}