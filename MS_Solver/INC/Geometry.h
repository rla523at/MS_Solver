#pragma once
#include "Reference_Geometry.h"

class Geometry
{
public:
	Geometry(const std::shared_ptr<const Reference_Geometry>& reference_goemetry, std::vector<Euclidean_Vector>&& consisting_nodes);

public://Command
	void change_points(std::vector<Euclidean_Vector>&& nodes);

public://Query
	bool operator==(const Geometry& other) const;

	Euclidean_Vector center_point(void) const { return this->mapping_vf_(this->reference_geometry_->center_point()); };
	std::vector<Geometry> face_geometries(void) const;
	const Quadrature_Rule& get_quadrature_rule(const ushort integrand_order) const;
	bool is_line(void) const;
	bool is_simplex(void) const;
	ushort num_post_nodes(const ushort post_order) const;
	ushort num_post_elements(const ushort post_order) const;
	ushort num_vertices(void) const;
	Euclidean_Vector normalized_normal_vector(const Euclidean_Vector& node) const;
	std::vector<Euclidean_Vector> normalized_normal_vectors(const std::vector<Euclidean_Vector>& points) const;
	Vector_Function<Polynomial> orthonormal_basis_vector_function(const ushort solution_order) const;
	std::vector<Euclidean_Vector> post_points(const ushort post_order) const;
	std::vector<Euclidean_Vector> post_element_centers(const ushort post_order) const;
	std::vector<std::vector<int>> post_connectivities(const ushort post_order, const size_t connectivity_start_index) const;
	std::vector<double> projected_volumes(void) const;
	Euclidean_Vector random_point(void) const { return this->mapping_vf_(this->reference_geometry_->random_point()); }
	std::vector<std::vector<Euclidean_Vector>> set_of_face_points(void) const;
	std::vector<std::vector<Euclidean_Vector>> set_of_sub_simplex_vertices(void) const;
	std::vector<Geometry> sub_simplex_geometries(void) const;
	double volume(void) const;
	std::vector<Euclidean_Vector> vertices(void) const;

protected:	
	//ushort check_space_dimension(void) const;
	Vector_Function<Polynomial> initial_basis_vector_function(const ushort solution_order) const;
	bool is_axis_parallel_node(const Euclidean_Vector& node) const;
	bool is_on_axis_plane(const ushort axis_tag) const;
	bool is_on_this_axis_plane(const ushort axis_tag, const double reference_value) const;
	bool is_on_same_axis_plane(const Geometry& other) const;
	//Vector_Function<Polynomial> make_mapping_function(void) const;
	Quadrature_Rule make_quadrature_rule(const ushort integrand_order) const;

protected:
	ushort space_dimension_ = 0;
	std::shared_ptr<const Reference_Geometry> reference_geometry_;
	std::vector<Euclidean_Vector> points_;
	Vector_Function<Polynomial> mapping_vf_;
	Vector_Function<Polynomial> normal_vf_;
	Irrational_Function scale_f_;
	mutable std::map<size_t, Quadrature_Rule> degree_to_quadrature_rule_;
};

#include <concepts>
namespace ms
{	
	template <typename T, typename Container>	std::vector<T> extract_by_index(const std::vector<T>& set, const Container& indexes) 
	{
		const auto num_extracted_value = indexes.size();

		std::vector<T> extracted_values;
		extracted_values.reserve(num_extracted_value);

		for (const auto& index : indexes)
			extracted_values.push_back(set[index]);

		return extracted_values;
	}
	template <typename F> double integrate(const F& integrand, const Quadrature_Rule& quadrature_rule)
	{
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += integrand(QP_set[i]) * QW_set[i];

		return result;
	}
	template <typename F> double integrate(const F& integrand, const Geometry& geometry)
	{
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand.degree());
		return ms::integrate(integrand, quadrature_rule);
	}
	double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry);	
	double L2_Norm(const Polynomial& function, const Geometry& geometry);
	Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& functions, const Geometry& geometry);	
}