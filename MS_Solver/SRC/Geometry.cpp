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