#include "../INC/Geometry.h"

Geometry::Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes) {
	this->reference_geometry_ = Reference_Geometry_Factory::make(figure, order);
	this->nodes_ = std::move(consisting_nodes);
	this->space_dimension_ = this->check_space_dimension();
	this->mapping_function_ = this->make_mapping_function();
	this->scale_function_ = this->reference_geometry_->scale_function(mapping_function_);
}

Euclidean_Vector Geometry::center_node(void) const {
	return this->mapping_function_(this->reference_geometry_->center_node());
}

ushort Geometry::num_post_nodes(const ushort post_order) const {
	return this->reference_geometry_->num_post_nodes(post_order);
}

ushort Geometry::num_post_elements(const ushort post_order) const {
	return this->reference_geometry_->num_post_elements(post_order);
}

std::vector<Euclidean_Vector> Geometry::post_nodes(const ushort post_order) const {
	const auto& ref_post_nodes = this->reference_geometry_->get_post_nodes(post_order);
	const auto num_post_nodes = ref_post_nodes.size();

	std::vector<Euclidean_Vector> post_nodes;
	post_nodes.reserve(num_post_nodes);

	for (const auto& ref_post_node : ref_post_nodes)
		post_nodes.push_back(this->mapping_function_(ref_post_node));
	
	return post_nodes;
}

std::vector<std::vector<int>> Geometry::post_connectivities(const ushort post_order, const size_t connectivity_start_index) const {
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

Vector_Function<Polynomial> Geometry::orthonormal_basis_vector_function(const ushort solution_order) const {
	const auto initial_basis_vector_function = this->initial_basis_vector_function(solution_order);
	return ms::Gram_Schmidt_process(initial_basis_vector_function, *this);
}


double Geometry::volume(void) const {
	const auto& quadrature_rule = this->get_quadrature_rule(0);

	auto volume = 0.0;
	for (const auto weight : quadrature_rule.weights)
		volume += weight;

	return volume;
}

const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort integrand_order) const {
	if (this->integrand_order_to_quadrature_rule_.find(integrand_order) == this->integrand_order_to_quadrature_rule_.end())
		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->make_quadrature_rule(integrand_order));

	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
}

Vector_Function<Polynomial> Geometry::initial_basis_vector_function(const ushort solution_order) const {
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

ushort Geometry::check_space_dimension(void) const {
	const auto expect_dimension = static_cast<ushort>(this->nodes_.front().size());

	const auto num_nodes = this->nodes_.size();
	for (ushort i = 1; i < num_nodes; ++i)
		REQUIRE(expect_dimension == this->nodes_[i].size(), "every node should have same space dimension");

	return expect_dimension;
}

Vector_Function<Polynomial> Geometry::make_mapping_function(void) const {
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

Quadrature_Rule Geometry::make_quadrature_rule(const ushort integrand_order) const {
	const auto reference_integrand_order = integrand_order + this->reference_geometry_->scale_function_order();
	const auto& ref_quadrature_rule = this->reference_geometry_->get_quadrature_rule(reference_integrand_order);
	
	const auto num_QP = ref_quadrature_rule.nodes.size();
	std::vector<Euclidean_Vector> transformed_QP;
	transformed_QP.reserve(num_QP);

	for (const auto& ref_node : ref_quadrature_rule.nodes)
		transformed_QP.push_back(this->mapping_function_(ref_node));

	std::vector<double> transformed_QW(num_QP);

	for (size_t i = 0; i < num_QP; ++i) {
		const auto& point = ref_quadrature_rule.nodes[i];
		const auto& weight = ref_quadrature_rule.weights[i];

		transformed_QW[i] = this->scale_function_(point) * weight;
	}

	return { transformed_QP, transformed_QW };
}

namespace ms
{
	double integrate(const Polynomial& integrand, const Quadrature_Rule& quadrature_rule) {
		const auto& QP_set = quadrature_rule.nodes;
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
		const auto& QP_set = quadrature_rule.nodes;
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
		const auto range_dimension = functions.range_dimension();

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