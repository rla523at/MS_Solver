#include "../INC/Geometry.h"

Geometry::Geometry(std::unique_ptr<Reference_Geometry>&& reference_goemetry, std::vector<Euclidean_Vector>&& consisting_points)
	: reference_geometry_(std::move(reference_goemetry)),
	  points_(std::move(consisting_points))
{
	this->space_dimension_ = this->check_space_dimension();
	this->mapping_vf_ = this->make_mapping_function();
	this->normal_vf_ = this->reference_geometry_->normal_vector_function(this->mapping_vf_);
	this->scale_f_ = this->reference_geometry_->scale_function(this->mapping_vf_);
};

void Geometry::change_points(std::vector<Euclidean_Vector>&& points)
{
	this->points_ = std::move(points);
	this->mapping_vf_ = this->make_mapping_function();
	this->normal_vf_ = this->reference_geometry_->normal_vector_function(this->mapping_vf_);
	this->scale_f_ = this->reference_geometry_->scale_function(this->mapping_vf_);
}

bool Geometry::operator==(const Geometry& other) const
{
	return this->points_ == other.points_ 
		&& *this->reference_geometry_ == *other.reference_geometry_;
}

Euclidean_Vector Geometry::center_point(void) const
{
	return this->mapping_vf_(this->reference_geometry_->center_point());
}

std::vector<Geometry> Geometry::face_geometries(void) const
{
	auto set_of_face_points = this->set_of_face_points();
	auto face_reference_geometries = this->reference_geometry_->face_reference_geometries();
	const auto num_face = face_reference_geometries.size();

	std::vector<Geometry> face_geometries;
	face_geometries.reserve(num_face);

	for (ushort i = 0; i < num_face; ++i)
	{
		face_geometries.push_back({ std::move(face_reference_geometries[i]), std::move(set_of_face_points[i]) });
	}

	return face_geometries;
}

ushort Geometry::num_post_nodes(const ushort post_order) const
{
	return this->reference_geometry_->num_post_nodes(post_order);
}

ushort Geometry::num_post_elements(const ushort post_order) const
{
	return this->reference_geometry_->num_post_elements(post_order);
}

ushort Geometry::num_vertices(void) const
{
	return this->reference_geometry_->num_vertices();
}

Euclidean_Vector Geometry::normalized_normal_vector(const Euclidean_Vector& point) const
{
	auto normal_v = this->normal_vf_(point);
	normal_v.normalize();
	return normal_v;
}

std::vector<Euclidean_Vector> Geometry::normalized_normal_vectors(const std::vector<Euclidean_Vector>& points) const
{
	const auto num_points = points.size();

	std::vector<Euclidean_Vector> normalized_normal_vectors(num_points);
	for (ushort i = 0; i < num_points; ++i)
	{
		normalized_normal_vectors[i] = this->normalized_normal_vector(points[i]);
	}

	return normalized_normal_vectors;
}

std::vector<Euclidean_Vector> Geometry::post_points(const ushort post_order) const
{
	const auto& ref_post_points = this->reference_geometry_->get_post_points(post_order);
	const auto num_post_points = ref_post_points.size();

	std::vector<Euclidean_Vector> post_points;
	post_points.reserve(num_post_points);

	for (const auto& point : ref_post_points)
	{
		post_points.push_back(this->mapping_vf_(point));
	}

	return post_points;
}

std::vector<Euclidean_Vector> Geometry::post_element_centers(const ushort post_order) const
{
	const auto post_points = this->post_points(post_order);
	const auto& ref_connectivities = this->reference_geometry_->get_connectivities(post_order);

	const auto num_post_element = ref_connectivities.size();
	std::vector<Euclidean_Vector> post_element_centers;
	post_element_centers.reserve(num_post_element);

	for (const auto& connectivity : ref_connectivities)
	{
		Euclidean_Vector center(this->space_dimension_);

		for (const auto index : connectivity)
		{
			center += post_points[index];
		}

		center *= 1.0 / connectivity.size();

		post_element_centers.push_back(center);
	}

	return post_element_centers;
}

std::vector<std::vector<int>> Geometry::post_connectivities(const ushort post_order, const size_t connectivity_start_index) const
{
	const auto& ref_connectivities = this->reference_geometry_->get_connectivities(post_order);

	const auto num_connectivity = ref_connectivities.size();
	std::vector<std::vector<int>> connectivities(num_connectivity);

	for (ushort i = 0; i < num_connectivity; ++i)
	{
		auto& connectivity = connectivities[i];
		const auto& ref_connecitivity = ref_connectivities[i];

		const auto num_index = ref_connecitivity.size();
		connectivity.resize(num_index);

		for (ushort j = 0; j < num_index; ++j)
		{
			const auto new_index = static_cast<int>(ref_connecitivity[j] + connectivity_start_index);
			connectivity[j] = new_index;
		}
	}

	return connectivities;
}

Vector_Function<Polynomial> Geometry::orthonormal_basis_vector_function(const ushort solution_order) const
{
	const auto initial_basis_vector_function = this->initial_basis_vector_function(solution_order);
	return ms::Gram_Schmidt_process(initial_basis_vector_function, *this);
}

std::vector<double> Geometry::projected_volumes(void) const
{
	//This only work for linear mesh and convex geometry and cell

	if (this->space_dimension_ == 2)
	{
		double x_projected_volume = 0.0;
		double y_projected_volume = 0.0;

		const auto set_of_face_points = this->set_of_face_points();
		for (const auto& points : set_of_face_points)
		{
			const auto& start_point_v = points[0];
			const auto& end_point_v = points[1];
			const auto distance_v = end_point_v - start_point_v;

			x_projected_volume += std::abs(distance_v.at(0));
			y_projected_volume += std::abs(distance_v.at(1));
		}

		return { 0.5 * y_projected_volume, 0.5 * x_projected_volume };
	}
	else if (this->space_dimension_ == 3)
	{
		double yz_projected_volume = 0.0;
		double xz_projected_volume = 0.0;
		double xy_projected_volume = 0.0;

		const auto face_geometries = this->face_geometries();
		for (const auto& geometry : face_geometries)
		{
			const auto normal_v = geometry.normalized_normal_vector(geometry.center_point());

			Euclidean_Vector yz_plane_normalized_normal_vector = { 1,0,0 };
			Euclidean_Vector xz_plane_normalized_normal_vector = { 0,1,0 };
			Euclidean_Vector xy_plane_normalized_normal_vector = { 0,0,1 };

			const auto volume = geometry.volume();

			yz_projected_volume += volume * std::abs(normal_v.inner_product(yz_plane_normalized_normal_vector));
			xz_projected_volume += volume * std::abs(normal_v.inner_product(xz_plane_normalized_normal_vector));
			xy_projected_volume += volume * std::abs(normal_v.inner_product(xy_plane_normalized_normal_vector));
		}

		return { 0.5 * yz_projected_volume, 0.5 * xz_projected_volume, 0.5 * xy_projected_volume };
	}
	else
	{
		EXCEPTION("not supproted space dimension");
		return {};
	}
}

std::vector<std::vector<Euclidean_Vector>> Geometry::set_of_face_points(void) const
{
	const auto set_of_face_point_index_sequences = this->reference_geometry_->set_of_face_node_index_sequences();
	const auto num_face = set_of_face_point_index_sequences.size();

	std::vector<std::vector<Euclidean_Vector>> set_of_face_points;
	set_of_face_points.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i)
	{
		set_of_face_points.push_back(ms::extract_by_index(this->points_, set_of_face_point_index_sequences[i]));
	}

	return set_of_face_points;
}

std::vector<std::vector<Euclidean_Vector>> Geometry::set_of_sub_simplex_vertices(void) const
{
	const auto set_of_sub_simplex_vertex_index_sequences = this->reference_geometry_->set_of_sub_simplex_vertex_index_sequences();
	const auto num_sub_simplex = set_of_sub_simplex_vertex_index_sequences.size();

	std::vector<std::vector<Euclidean_Vector>> set_of_simplex_vertices;
	set_of_simplex_vertices.reserve(num_sub_simplex);

	for (size_t i = 0; i < num_sub_simplex; ++i)
	{
		set_of_simplex_vertices.push_back(ms::extract_by_index(this->points_, set_of_sub_simplex_vertex_index_sequences[i]));
	}

	return set_of_simplex_vertices;
}

std::vector<Geometry> Geometry::sub_simplex_geometries(void) const
{
	auto sub_simplex_reference_geometries = this->reference_geometry_->sub_simplex_reference_geometries();
	auto set_of_sub_simplex_vertices = this->set_of_sub_simplex_vertices();
	
	const auto num_sub_simplex = sub_simplex_reference_geometries.size();

	std::vector<Geometry> sub_simplex_geometries;
	sub_simplex_geometries.reserve(num_sub_simplex);

	for (ushort i = 0; i < num_sub_simplex; ++i)
	{
		sub_simplex_geometries.push_back({ std::move(sub_simplex_reference_geometries[i]), std::move(set_of_sub_simplex_vertices[i]) });
	}

	return sub_simplex_geometries;
}

double Geometry::volume(void) const
{
	const auto& quadrature_rule = this->get_quadrature_rule(0);

	auto volume = 0.0;
	for (const auto weight : quadrature_rule.weights)
	{
		volume += weight;
	}

	return volume;
}

std::vector<Euclidean_Vector> Geometry::vertices(void) const
{
	const auto num_vertices = this->reference_geometry_->num_vertices();
	return { this->points_.begin(), this->points_.begin() + num_vertices };
}

bool Geometry::is_line(void) const
{
	return this->reference_geometry_->is_line();
}

bool Geometry::is_simplex(void) const
{
	return this->reference_geometry_->is_simplex();
}

const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort integrand_degree) const
{
	if (!this->degree_to_quadrature_rule_.contains(integrand_degree))
		this->degree_to_quadrature_rule_.emplace(integrand_degree, this->make_quadrature_rule(integrand_degree));

	return this->degree_to_quadrature_rule_.at(integrand_degree);
}

ushort Geometry::check_space_dimension(void) const
{
	const auto expect_dimension = static_cast<ushort>(this->points_.front().size());

	const auto num_nodes = this->points_.size();
	for (ushort i = 1; i < num_nodes; ++i)
	{
		REQUIRE(expect_dimension == this->points_[i].size(), "every point should have same space dimension");
	}

	return expect_dimension;
}

Vector_Function<Polynomial> Geometry::initial_basis_vector_function(const ushort solution_order) const
{
	const auto num_basis = ms::combination_with_repetition(1 + this->space_dimension_, solution_order);

	std::vector<Polynomial> initial_basis_functions(num_basis);

	ushort index = 0;
	if (this->space_dimension_ == 2)
	{
		Polynomial x("x0");
		Polynomial y("x1");

		const auto center_point = this->center_point();
		const auto x_c = center_point.at(0);
		const auto y_c = center_point.at(1);

		//1 (x - x_c) (y - y_c)  ...
		for (ushort a = 0; a <= solution_order; ++a)
		{
			for (ushort b = 0; b <= a; ++b)
			{
				initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
			}
		}
	}
	else if (this->space_dimension_ == 3)
	{
		Polynomial x("x0");
		Polynomial y("x1");
		Polynomial z("x2");

		const auto center_point = this->center_point();
		const auto x_c = center_point.at(0);
		const auto y_c = center_point.at(1);
		const auto z_c = center_point.at(2);

		//1 (x - x_c) (y - y_c) (z - z_c) ...
		for (ushort a = 0; a <= solution_order; ++a)
		{
			for (ushort b = 0; b <= a; ++b)
			{
				for (ushort c = 0; c <= b; ++c)
				{
					initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ (b - c)) * ((z - z_c) ^ c);
				}
			}
		}
	}
	else
		EXCEPTION("not supported space dimension");

	return initial_basis_functions;
}

bool Geometry::is_axis_parallel_node(const Euclidean_Vector& node) const
{
	for (const auto& my_node : this->points_)
	{
		if (my_node.is_axis_translation(node))
		{
			return true;
		}
	}
	return false;
}

bool Geometry::is_on_axis_plane(const ushort axis_tag) const
{
	if (this->space_dimension_ <= axis_tag)
	{
		return false;
	}

	const auto& ref_node = this->points_.front();

	constexpr auto epsilon = 1.0E-10;
	for (ushort j = 1; j < this->points_.size(); ++j)
	{
		if (std::abs(this->points_[j][axis_tag] - ref_node[axis_tag]) > epsilon)
		{
			return false;
		}
	}
	return true;
}

bool Geometry::is_on_this_axis_plane(const ushort axis_tag, const double reference_value) const
{
	if (this->space_dimension_ <= axis_tag)
	{
		return false;
	}

	constexpr auto epsilon = 1.0E-10;
	for (const auto& node : this->points_)
	{
		if (std::abs(node[axis_tag] - reference_value) > epsilon)
		{
			return false;
		}
	}

	return true;
}


bool Geometry::is_on_same_axis_plane(const Geometry& other) const
{
	for (ushort i = 0; i < this->space_dimension_; ++i)
	{
		if (this->is_on_axis_plane(i))
		{
			const auto reference_value = this->points_.front()[i];

			if (other.is_on_this_axis_plane(i, reference_value))
			{
				return true;
			}
		}
	}

	return false;
}

Vector_Function<Polynomial> Geometry::make_mapping_function(void) const
{
	//	X = CM
	//	X : mapped node matrix			
	//	C : mapping coefficient matrix	
	//	M : mapping monomial matrix

	const auto num_nodes = this->points_.size();
	Matrix X(this->space_dimension_, num_nodes);
	for (size_t j = 0; j < num_nodes; ++j)
	{
		X.change_column(j, this->points_[j]);
	}

	const auto& inv_M = this->reference_geometry_->get_inverse_mapping_monomial_matrix();

	const auto C = X * inv_M;

	const auto& monomial_vf = this->reference_geometry_->get_mapping_monomial_vector_function();

	return C * monomial_vf;
}

Quadrature_Rule Geometry::make_quadrature_rule(const ushort integrand_order) const
{
	const auto reference_integrand_order = integrand_order + this->reference_geometry_->scale_function_order();
	const auto& ref_quadrature_rule = this->reference_geometry_->get_quadrature_rule(reference_integrand_order);

	const auto num_QP = ref_quadrature_rule.points.size();
	std::vector<Euclidean_Vector> transformed_QPs;
	transformed_QPs.reserve(num_QP);

	for (const auto& ref_point_v : ref_quadrature_rule.points)
	{
		transformed_QPs.push_back(this->mapping_vf_(ref_point_v));
	}

	std::vector<double> transformed_QW(num_QP);

	for (size_t i = 0; i < num_QP; ++i)
	{
		const auto& point_v = ref_quadrature_rule.points[i];
		const auto& weight = ref_quadrature_rule.weights[i];

		transformed_QW[i] = this->scale_f_(point_v) * weight;
	}

	return { transformed_QPs, transformed_QW };
}

namespace ms
{
	double integrate(const Polynomial& integrand, const Quadrature_Rule& quadrature_rule)
	{
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += integrand(QP_set[i]) * QW_set[i];

		return result;
	}

	double integrate(const Polynomial& integrand, const Geometry& geometry)
	{
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand.degree());
		return ms::integrate(integrand, quadrature_rule);
	}

	double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry)
	{
		const auto integrand_degree = f1.degree() + f2.degree();

		const auto quadrature_rule = geometry.get_quadrature_rule(integrand_degree);
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += f1(QP_set[i]) * f2(QP_set[i]) * QW_set[i];

		return result;
	}

	double L2_Norm(const Polynomial& function, const Geometry& geometry)
	{
		return std::sqrt(ms::inner_product(function, function, geometry));
	}

	Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& vf, const Geometry& geometry)
	{
		const auto range_dimension = vf.size();

		std::vector<Polynomial> normalized_functions(range_dimension);

		for (ushort i = 0; i < range_dimension; ++i)
		{
			normalized_functions[i] = vf[i];

			for (ushort j = 0; j < i; ++j)
				normalized_functions[i] -= ms::inner_product(normalized_functions[i], normalized_functions[j], geometry) * normalized_functions[j];

			normalized_functions[i] *= 1.0 / ms::L2_Norm(normalized_functions[i], geometry);
		}

		return normalized_functions;
	}
}