#include "../INC/Geometry.h"

Geometry::Geometry(const Figure figure, const ushort order, std::vector<Euclidean_Vector>&& consisting_nodes) {
	this->reference_geometry_ = Reference_Geometry_Factory::make(figure, order);
	this->nodes_ = std::move(consisting_nodes);
	this->space_dimension_ = this->check_space_dimension();
	this->mapping_function_ = this->make_mapping_function();
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

const Quadrature_Rule& Geometry::get_quadrature_rule(const ushort integrand_order) const {
	if (this->integrand_order_to_quadrature_rule_.find(integrand_order) == this->integrand_order_to_quadrature_rule_.end())
		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->make_quadrature_rule(integrand_order));

	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
}

Vector_Function<Polynomial> Geometry::orthonormal_basis_vector_function(const ushort solution_order) const {
	const auto initial_basis_set = this->initial_basis_vector_function(solution_order);
	return ms::Gram_Schmidt_process(initial_basis_set, *this);
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
	const auto scale_function = this->scale_function(mapping_function);
	const auto reference_integrand_order = integrand_order + this->scale_function_order();

	const auto& ref_quadrature_rule = this->reference_geometry_->get_quadrature_rule(reference_integrand_order);
	
	const auto transformed_QP = mapping_function(ref_quadrature_rule.nodes);

	const auto num_QP = transformed_QP.size();
	std::vector<double> transformed_QW(num_QP);

	for (size_t i = 0; i < num_QP; ++i) {
		const auto& point = ref_quadrature_rule.nodes[i];
		const auto& weight = ref_quadrature_rule.weights[i];

		transformed_QW[i] = scale_function(point) * weight;
	}

	return { transformed_QP, transformed_QW };
}

Irrational_Function Geometry::make_scale_function(const Vector_Function<Polynomial>& mapping_function) const {
	if constexpr (space_dimension == 2) {
		switch (this->figure_) {
		case Figure::line: {
			constexpr ushort r = 0;
			const auto mf_r = mapping_function.differentiate(r);
			return mf_r.L2_norm();
		}
		case Figure::triangle:
		case Figure::quadrilateral: {
			constexpr ushort r = 0;
			constexpr ushort s = 1;
			const auto mf_r = mapping_function.differentiate(r);
			const auto mf_s = mapping_function.differentiate(s);
			const auto cross_product = mf_r.cross_product(mf_s);
			return cross_product.L2_norm();
		}
		default:
			throw std::runtime_error("not supported figure");
			return Irrational_Function();
		}
	}
	else if constexpr (space_dimension == 3) {
		switch (this->figure_) {
		case Figure::triangle:
		case Figure::quadrilateral: {
			constexpr ushort r = 0;
			constexpr ushort s = 1;
			const auto mf_r = mapping_function.differentiate(r);
			const auto mf_s = mapping_function.differentiate(s);
			const auto cross_product = mf_r.cross_product(mf_s);
			return cross_product.L2_norm();
		}
		case Figure::tetrahedral:
		case Figure::hexahedral: {
			constexpr ushort r = 0;
			constexpr ushort s = 1;
			constexpr ushort t = 2;
			const auto mf_r = mapping_function.differentiate(r);
			const auto mf_s = mapping_function.differentiate(s);
			const auto mf_t = mapping_function.differentiate(t);

			return ms::scalar_triple_product(mf_r, mf_s, mf_t).be_absolute();
		}
		default:
			throw std::runtime_error("not supported figure");
			return Irrational_Function();
		}
	}
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