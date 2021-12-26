#include "../INC/Reference_Geometry.h"

bool Reference_Geometry::operator==(const Reference_Geometry& other) const
{
	return (this->figure() == other.figure()) && (this->order_ == other.order_);
}

bool Reference_Geometry::operator!=(const Reference_Geometry& other) const
{
	return !(*this == other);
}

const std::vector<Euclidean_Vector>& Reference_Geometry::get_mapping_nodes(void) const
{
	const auto figure_index = static_cast<int>(this->figure());
	return this->figure_index_and_order_to_mapping_points_table_[figure_index][this->order_];
}

const Vector_Function<Polynomial>& Reference_Geometry::get_mapping_monomial_vector_function(void) const
{
	const auto figure_index = static_cast<int>(this->figure());
	return this->figure_index_and_order_to_mapping_monomial_vf_[figure_index][this->order_];
}

const Matrix& Reference_Geometry::get_inverse_mapping_monomial_matrix(void) const
{
	const auto figure_index = static_cast<int>(this->figure());
	return this->figure_index_and_order_to_inverse_mapping_monomial_m_[figure_index][this->order_];
}

const std::vector<Euclidean_Vector>& Reference_Geometry::get_post_points(const ushort post_order) const
{
	const auto figure_index = static_cast<int>(this->figure());

	if (this->set_of_post_nodes_[figure_index].size() <= post_order)
	{
		const auto new_size = post_order + 1;
		this->set_of_post_nodes_[figure_index].resize(new_size);
	}

	if (this->set_of_post_nodes_[figure_index][post_order].empty())
	{
		this->set_of_post_nodes_[figure_index][post_order] = this->make_post_nodes(post_order);
	}

	return this->set_of_post_nodes_[figure_index][post_order];
}

const std::vector<std::vector<uint>>& Reference_Geometry::get_connectivities(const ushort post_order) const
{
	const auto figure_index = static_cast<int>(this->figure());

	if (this->set_of_connectivities_[figure_index].size() <= post_order)
	{
		const auto new_size = post_order + 1;
		this->set_of_connectivities_[figure_index].resize(new_size);
	}

	if (this->set_of_connectivities_[figure_index][post_order].empty())
	{
		this->set_of_connectivities_[figure_index][post_order] = this->make_connectivities(post_order);
	}

	return this->set_of_connectivities_[figure_index][post_order];
}

const Quadrature_Rule& Reference_Geometry::get_quadrature_rule(const ushort integrand_order) const
{
	const auto figure_index = static_cast<int>(this->figure());

	if (this->quadrature_rules_[figure_index].size() <= integrand_order)
	{
		const auto new_size = integrand_order + 1;
		this->quadrature_rules_[figure_index].resize(new_size);
	}

	if (this->quadrature_rules_[figure_index][integrand_order].is_empty())
	{
		this->quadrature_rules_[figure_index][integrand_order] = this->make_quadrature_rule(integrand_order);
	}

	return this->quadrature_rules_[figure_index][integrand_order];
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
	const auto figure_index = static_cast<int>(this->figure());

	if (this->figure_index_and_order_to_mapping_points_table_[figure_index].size() < this->order_)
	{
		const auto new_size = this->order_ + 1;
		this->figure_index_and_order_to_mapping_points_table_[figure_index].resize(new_size);
		this->figure_index_and_order_to_mapping_monomial_vf_[figure_index].resize(new_size);
		this->figure_index_and_order_to_inverse_mapping_monomial_m_[figure_index].resize(new_size);
	}

	if (this->figure_index_and_order_to_mapping_points_table_[figure_index][this->order_].empty())
	{
		this->figure_index_and_order_to_mapping_points_table_[figure_index][this->order_] = this->make_mapping_nodes();
		this->figure_index_and_order_to_mapping_monomial_vf_[figure_index][this->order_] = this->make_mapping_monomial_vector_function();
		this->figure_index_and_order_to_inverse_mapping_monomial_m_[figure_index][this->order_] = this->make_inverse_mapping_monomial_matrix();
	}
}

Matrix Reference_Geometry::make_inverse_mapping_monomial_matrix(void) const
{
	const auto& mapping_points = this->get_mapping_nodes();
	const auto& mapping_monomial_vf = this->get_mapping_monomial_vector_function();

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