#include "../INC/Geometry.h"

Geometry::Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes) 
{
	this->reference_geometry_ = Reference_Geometry_Factory::make(figure, order);
	this->nodes_ = std::move(consisting_nodes);
	this->space_dimension_ = this->check_space_dimension();
	this->mapping_function_ = this->make_mapping_function();
	this->scale_function_ = this->reference_geometry_->scale_function(mapping_function_);
}

Geometry::Geometry(std::unique_ptr<Reference_Geometry>&& reference_goemetry, std::vector<Euclidean_Vector>&& consisting_nodes)
	:	reference_geometry_(std::move(reference_goemetry)),
		nodes_(std::move(consisting_nodes)) 
{
	this->space_dimension_ = this->check_space_dimension();
	this->mapping_function_ = this->make_mapping_function();
	this->scale_function_ = this->reference_geometry_->scale_function(mapping_function_);
};

void Geometry::change_nodes(std::vector<Euclidean_Vector>&& nodes)
{
	*this = Geometry(std::move(this->reference_geometry_), std::move(nodes));
}

bool Geometry::operator==(const Geometry& other) const
{
	return this->nodes_ == other.nodes_ &&
		*this->reference_geometry_ == *other.reference_geometry_;
}

Euclidean_Vector Geometry::center_node(void) const 
{
	return this->mapping_function_(this->reference_geometry_->center_node());
}

std::vector<Geometry> Geometry::face_geometries(void) const
{
	auto set_of_face_nodes = this->set_of_face_nodes();	
	auto face_reference_geometries = this->reference_geometry_->face_reference_geometries();
	const auto num_face = face_reference_geometries.size();
	
	std::vector<Geometry> face_geometries;
	face_geometries.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i) 
	{
		face_geometries.push_back({ std::move(face_reference_geometries[i]), std::move(set_of_face_nodes[i]) });
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

Euclidean_Vector Geometry::normalized_normal_vector(const Euclidean_Vector& node) const 
{
	const auto normal_vector_function = this->reference_geometry_->normal_vector_function(this->mapping_function_);
	auto normal_vector = Euclidean_Vector(normal_vector_function(node));

	return normal_vector.normalize();
}

std::vector<Euclidean_Vector> Geometry::normalized_normal_vectors(const std::vector<Euclidean_Vector>& points) const
{
	const auto normal_vector_function = this->reference_geometry_->normal_vector_function(this->mapping_function_);

	const auto num_points = points.size();

	std::vector<Euclidean_Vector> normalized_normal_vectors;
	normalized_normal_vectors.reserve(num_points);

	for (int i = 0; i < num_points; ++i)
	{
		auto normal_vector = Euclidean_Vector(normal_vector_function(points[i]));
		normalized_normal_vectors[i] = normal_vector.normalize();
	}

	return normalized_normal_vectors;
}


std::vector<Euclidean_Vector> Geometry::post_nodes(const ushort post_order) const 
{
	const auto& ref_post_nodes = this->reference_geometry_->get_post_nodes(post_order);
	const auto num_post_nodes = ref_post_nodes.size();

	std::vector<Euclidean_Vector> post_nodes;
	post_nodes.reserve(num_post_nodes);

	for (const auto& ref_post_node : ref_post_nodes)
		post_nodes.push_back(this->mapping_function_(ref_post_node));
	
	return post_nodes;
}

std::vector<std::vector<int>> Geometry::post_connectivities(const ushort post_order, const size_t connectivity_start_index) const 
{
	const auto& ref_connectivities = this->reference_geometry_->get_connectivities(post_order);

	const auto num_connectivity = ref_connectivities.size();
	std::vector<std::vector<int>> connectivities(num_connectivity);

	for (ushort i = 0; i < num_connectivity; ++i) {
		auto& connectivity = connectivities[i];
		const auto& ref_connecitivity = ref_connectivities[i];

		const auto num_index = ref_connecitivity.size();
		connectivity.resize(num_index);

		for (ushort j = 0; j < num_index; ++j) {
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

std::vector<double> Geometry::projected_volume(void) const 
{
	//This only work for linear mesh and convex geometry.

	if (this->space_dimension_ == 2)
	{
		double x_projected_volume = 0.0;
		double y_projected_volume = 0.0;

		const auto set_of_face_nodes = this->set_of_face_nodes();
		for (const auto& face_nodes : set_of_face_nodes) {
			const auto& start_node = face_nodes[0];
			const auto& end_node = face_nodes[1];
			const auto node_to_node = end_node - start_node;

			x_projected_volume += std::abs(node_to_node.at(0));
			y_projected_volume += std::abs(node_to_node.at(1));
		}

		return { 0.5 * y_projected_volume, 0.5 * x_projected_volume };
	}
	else if(this->space_dimension_ == 3) 
	{
		double yz_projected_volume = 0.0;
		double xz_projected_volume = 0.0;
		double xy_projected_volume = 0.0;

		const auto face_geometries = this->face_geometries();
		for (const auto& geometry : face_geometries) {
			const auto normal_vector = geometry.normalized_normal_vector(geometry.center_node());

			Euclidean_Vector yz_plane_normalized_normal_vector = { 1,0,0 };
			Euclidean_Vector xz_plane_normalized_normal_vector = { 0,1,0 };
			Euclidean_Vector xy_plane_normalized_normal_vector = { 0,0,1 };

			const auto volume = geometry.volume();

			yz_projected_volume += volume * std::abs(normal_vector.inner_product(yz_plane_normalized_normal_vector));
			xz_projected_volume += volume * std::abs(normal_vector.inner_product(xz_plane_normalized_normal_vector));
			xy_projected_volume += volume * std::abs(normal_vector.inner_product(xy_plane_normalized_normal_vector));
		}

		return { 0.5 * yz_projected_volume, 0.5 * xz_projected_volume, 0.5 * xy_projected_volume };
	}
	else 
	{
		EXCEPTION("not supproted space dimension");
		return {};
	}
}

std::vector<std::vector<Euclidean_Vector>> Geometry::set_of_face_nodes(void) const 
{
	const auto set_of_face_node_index_sequences = this->reference_geometry_->set_of_face_node_index_sequences();
	const auto num_face = set_of_face_node_index_sequences.size();

	std::vector<std::vector<Euclidean_Vector>> set_of_face_nodes;
	set_of_face_nodes.reserve(num_face);

	for (size_t i = 0; i < num_face; ++i)
		set_of_face_nodes.push_back(ms::extract_by_index(this->nodes_, set_of_face_node_index_sequences[i]));

	return set_of_face_nodes;
}

double Geometry::volume(void) const 
{
	const auto& quadrature_rule = this->get_quadrature_rule(0);

	auto volume = 0.0;
	for (const auto weight : quadrature_rule.weights)
		volume += weight;

	return volume;
}

bool Geometry::is_line(void) const
{
	return this->reference_geometry_->is_line();
}

const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort integrand_order) const 
{
	if (this->integrand_order_to_quadrature_rule_.find(integrand_order) == this->integrand_order_to_quadrature_rule_.end())
		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->make_quadrature_rule(integrand_order));

	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
}

ushort Geometry::check_space_dimension(void) const 
{
	const auto expect_dimension = static_cast<ushort>(this->nodes_.front().size());

	const auto num_nodes = this->nodes_.size();
	for (ushort i = 1; i < num_nodes; ++i)
		REQUIRE(expect_dimension == this->nodes_[i].size(), "every node should have same space dimension");

	return expect_dimension;
}

Vector_Function<Polynomial> Geometry::initial_basis_vector_function(const ushort solution_order) const
{
	const auto num_basis = ms::combination_with_repetition(1 + this->space_dimension_, solution_order);

	std::vector<Polynomial> initial_basis_functions(num_basis);

	ushort index = 0;
	if (this->space_dimension_ == 2) {
		Polynomial x("x0");
		Polynomial y("x1");

		const auto center_node = this->center_node();
		const auto x_c = center_node.at(0);
		const auto y_c = center_node.at(1);

		//1 (x - x_c) (y - y_c)  ...
		for (ushort a = 0; a <= solution_order; ++a)
			for (ushort b = 0; b <= a; ++b)
				initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ b);
	}
	else if (this->space_dimension_ == 3) {
		Polynomial x("x0");
		Polynomial y("x1");
		Polynomial z("x2");

		const auto center_node = this->center_node();
		const auto x_c = center_node.at(0);
		const auto y_c = center_node.at(1);
		const auto z_c = center_node.at(2);

		//1 (x - x_c) (y - y_c) (z - z_c) ...
		for (ushort a = 0; a <= solution_order; ++a)
			for (ushort b = 0; b <= a; ++b)
				for (ushort c = 0; c <= b; ++c)
					initial_basis_functions[index++] = ((x - x_c) ^ (a - b)) * ((y - y_c) ^ (b - c)) * ((z - z_c) ^ c);
	}
	else
		throw std::runtime_error("not supported space dimension");


	return initial_basis_functions;
}

bool Geometry::is_axis_parallel_node(const Euclidean_Vector& node) const 
{
	for (const auto& my_node : this->nodes_) 
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

	const auto& ref_node = this->nodes_.front();

	constexpr auto epsilon = 1.0E-10;
	for (ushort j = 1; j < this->nodes_.size(); ++j)
	{
		if (std::abs(this->nodes_[j][axis_tag] - ref_node[axis_tag]) > epsilon)
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
	for (const auto& node : this->nodes_)
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
			const auto reference_value = this->nodes_.front()[i];
			
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

	const auto num_nodes = this->nodes_.size();
	Matrix X(this->space_dimension_, num_nodes);
	for (size_t j = 0; j < num_nodes; ++j)
		X.change_column(j, this->nodes_[j]);

	const auto& inv_M = this->reference_geometry_->get_inverse_mapping_monomial_matrix();

	const auto C = X * inv_M;

	const auto& monomial_vector_function = this->reference_geometry_->get_mapping_monomial_vector_function();

	return C * monomial_vector_function;
}

Quadrature_Rule Geometry::make_quadrature_rule(const ushort integrand_order) const 
{
	const auto reference_integrand_order = integrand_order + this->reference_geometry_->scale_function_order();
	const auto& ref_quadrature_rule = this->reference_geometry_->get_quadrature_rule(reference_integrand_order);
	
	const auto num_QP = ref_quadrature_rule.points.size();
	std::vector<Euclidean_Vector> transformed_QP;
	transformed_QP.reserve(num_QP);

	for (const auto& ref_node : ref_quadrature_rule.points)
		transformed_QP.push_back(this->mapping_function_(ref_node));

	std::vector<double> transformed_QW(num_QP);

	for (size_t i = 0; i < num_QP; ++i) {
		const auto& point = ref_quadrature_rule.points[i];
		const auto& weight = ref_quadrature_rule.weights[i];

		transformed_QW[i] = this->scale_function_(point) * weight;
	}

	return { transformed_QP, transformed_QW };
}

namespace ms
{
	double integrate(const Polynomial& integrand, const Quadrature_Rule& quadrature_rule) {
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += integrand(QP_set[i]) * QW_set[i];

		return result;
	}

	double integrate(const Polynomial& integrand, const Geometry& geometry) {
		const auto quadrature_rule = geometry.get_quadrature_rule(integrand.degree());
		return ms::integrate(integrand, quadrature_rule);
	}

	double inner_product(const Polynomial& f1, const Polynomial& f2, const Geometry& geometry) {
		const auto integrand_degree = f1.degree() + f2.degree();

		const auto quadrature_rule = geometry.get_quadrature_rule(integrand_degree);
		const auto& QP_set = quadrature_rule.points;
		const auto& QW_set = quadrature_rule.weights;

		double result = 0.0;
		for (ushort i = 0; i < QP_set.size(); ++i)
			result += f1(QP_set[i]) * f2(QP_set[i]) * QW_set[i];

		return result;
	}

	double L2_Norm(const Polynomial& function, const Geometry& geometry) {
		return std::sqrt(ms::inner_product(function, function, geometry));
	}

	Vector_Function<Polynomial> Gram_Schmidt_process(const Vector_Function<Polynomial>& functions, const Geometry& geometry) {
		const auto range_dimension = functions.size();

		std::vector<Polynomial> normalized_functions(range_dimension);

		for (ushort i = 0; i < range_dimension; ++i) {
			normalized_functions[i] = functions[i];

			for (ushort j = 0; j < i; ++j)
				normalized_functions[i] -= ms::inner_product(normalized_functions[i], normalized_functions[j], geometry) * normalized_functions[j];

			normalized_functions[i] *= 1.0 / ms::L2_Norm(normalized_functions[i], geometry);
		}

		return normalized_functions;
	}
}