#include "../INC/Reference_Geometry.h"

ushort Reference_Geometry::check_space_dimension(const std::vector<Euclidean_Vector>& points) const
{
	const auto expect_dimension = static_cast<ushort>(points.front().size());

	const auto num_nodes = points.size();
	for (ushort i = 1; i < num_nodes; ++i)
	{
		REQUIRE(expect_dimension == points[i].size(), "every point should have same space dimension");
	}

	return expect_dimension;
}

bool Reference_Geometry::operator==(const Reference_Geometry& other) const
{
	return (this->figure() == other.figure()) && (this->order_ == other.order_);
}

bool Reference_Geometry::operator!=(const Reference_Geometry& other) const
{
	return !(*this == other);
}

const std::vector<Euclidean_Vector>& Reference_Geometry::get_post_points(const ushort post_order) const
{
	if (!this->post_order_to_post_points_.contains(post_order))
	{
		this->post_order_to_post_points_.emplace(post_order, this->make_post_points(post_order));
	}

	return this->post_order_to_post_points_.at(post_order);
}

const std::vector<std::vector<uint>>& Reference_Geometry::get_connectivities(const ushort post_order) const
{
	if (!this->post_order_to_connectivities_.contains(post_order))
	{
		this->post_order_to_connectivities_.emplace(post_order, this->make_connectivities(post_order));
	}

	return this->post_order_to_connectivities_.at(post_order);
}

const Quadrature_Rule& Reference_Geometry::get_quadrature_rule(const ushort integrand_order) const
{
	if (!this->integrand_order_to_quadrature_rule_.contains(integrand_order))
	{
		this->integrand_order_to_quadrature_rule_.emplace(integrand_order, this->make_quadrature_rule(integrand_order));
	}

	return this->integrand_order_to_quadrature_rule_.at(integrand_order);
}

Vector_Function<Polynomial> Reference_Geometry::make_mapping_function(const std::vector<Euclidean_Vector>& points) const
{
	//	X = CM
	//	X : mapped node matrix			
	//	C : mapping coefficient matrix	
	//	M : mapping monomial matrix

	const auto space_dimension = this->check_space_dimension(points);
	const auto num_points = points.size();
	Matrix X(space_dimension, num_points);
	for (size_t j = 0; j < num_points; ++j)
	{
		X.change_column(j, points[j]);
	}

	const auto& inv_M = this->order_to_inverse_mapping_monomial_m_.at(this->order_);

	const auto C = X * inv_M;

	const auto& monomial_vf = this->order_to_mapping_monomial_vf_.at(this->order_);

	return C * monomial_vf;
}

std::vector<ushort> Reference_Geometry::vertex_node_index_sequneces(void) const
{
	const auto num_vertices = this->num_vertices();
	std::vector<ushort> vnode_index_orders(num_vertices);

	for (ushort i = 0; i < num_vertices; ++i)
		vnode_index_orders[i] = i;

	return vnode_index_orders;
}

void Reference_Geometry::initialize(void)
{
	if (!this->order_to_mapping_points_.contains(this->order_))
	{
		this->order_to_mapping_points_.emplace(this->order_, this->make_mapping_points());
		this->order_to_mapping_monomial_vf_.emplace(this->order_, this->make_mapping_monomial_vector_function());
		this->order_to_inverse_mapping_monomial_m_.emplace(this->order_, this->make_inverse_mapping_monomial_matrix());
	}

	//const auto figure_index = static_cast<int>(this->figure());

	//if (this->figure_index_and_order_to_mapping_points_table_[figure_index].size() < this->order_)
	//{
	//	const auto new_size = this->order_ + 1;
	//	this->figure_index_and_order_to_mapping_points_table_[figure_index].resize(new_size);
	//	this->figure_index_and_order_to_mapping_monomial_vf_[figure_index].resize(new_size);
	//	this->figure_index_and_order_to_inverse_mapping_monomial_m_[figure_index].resize(new_size);
	//}

	//if (this->figure_index_and_order_to_mapping_points_table_[figure_index][this->order_].empty())
	//{
	//	this->figure_index_and_order_to_mapping_points_table_[figure_index][this->order_] = this->make_mapping_nodes();
	//	this->figure_index_and_order_to_mapping_monomial_vf_[figure_index][this->order_] = this->make_mapping_monomial_vector_function();
	//	this->figure_index_and_order_to_inverse_mapping_monomial_m_[figure_index][this->order_] = this->make_inverse_mapping_monomial_matrix();
	//}
}

Matrix Reference_Geometry::make_inverse_mapping_monomial_matrix(void) const
{
	const auto& mapping_points = this->order_to_mapping_points_.at(this->order_);
	const auto& mapping_monomial_vf = this->order_to_mapping_monomial_vf_.at(this->order_);

	//const auto& mapping_points = this->get_mapping_nodes();
	//const auto& mapping_monomial_vf = this->get_mapping_monomial_vector_function();

	const auto num_mapping_monomials = mapping_monomial_vf.size();
	Matrix mapping_monomial_matrix(num_mapping_monomials);

	for (size_t i = 0; i < num_mapping_monomials; ++i)
	{
		mapping_monomial_matrix.change_column(i, mapping_monomial_vf(mapping_points[i]));
	}

	mapping_monomial_matrix.inverse();

	return mapping_monomial_matrix;
}

std::vector<std::vector<uint>> Reference_Geometry::quadrilateral_connectivities(const std::array<uint, 4>& node_indexes) const
{
	//   2式式式式式3
	//   弛     弛   
	//   0式式式式式1
	return { { node_indexes[0],node_indexes[1],node_indexes[2] },
			 { node_indexes[1],node_indexes[3],node_indexes[2] } };
};