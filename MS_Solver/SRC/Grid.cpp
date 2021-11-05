#include "../INC/Grid.h"


size_t Grid::num_cells(void) const {
	return this->grid_elements_.cell_elements.size();
}

std::vector<Euclidean_Vector> Grid::cell_center_nodes(void) const {
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<Euclidean_Vector> center_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		center_nodes[i] = cell_elements[i].center_node();

	return center_nodes;
}

std::vector<ushort> Grid::cell_set_of_num_post_nodes(const ushort post_order) const {
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<ushort> cell_set_of_num_post_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		cell_set_of_num_post_nodes[i] = cell_elements[i].num_post_nodes(post_order);

	return cell_set_of_num_post_nodes;
}

std::vector<ushort> Grid::cell_set_of_num_post_elements(const ushort post_order) const{
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<ushort> cell_set_of_num_post_elements(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		cell_set_of_num_post_elements[i] = cell_elements[i].num_post_elements(post_order);

	return cell_set_of_num_post_elements;
}

std::vector<std::vector<Euclidean_Vector>> Grid::cell_set_of_post_nodes(const ushort post_order) const {
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		set_of_post_nodes[i] = cell_elements[i].post_nodes(post_order);

	return set_of_post_nodes;
}

std::vector<std::vector<int>> Grid::cell_set_of_connectivities(const ushort post_order, const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const {
	const auto& cell_elements = this->grid_elements_.cell_elements;

	size_t connectivity_start_index = 0;

	std::vector<std::vector<int>> set_of_connectivities;

	for (uint i = 0; i < cell_elements.size(); ++i) {
		auto icell_connectivities = cell_elements[i].post_connectivities(post_order, connectivity_start_index);
		ms::merge(set_of_connectivities, std::move(icell_connectivities));

		const auto& post_nodes = set_of_post_nodes[i];
		const auto num_post_nodes = post_nodes.size();
		connectivity_start_index += num_post_nodes;
	}

	return set_of_connectivities;
}