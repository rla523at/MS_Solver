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

	const auto num_cell = cell_elements.size();
	std::vector<std::vector<int>> set_of_connectivities;

	size_t connectivity_start_index = 0;

	for (uint i = 0; i < num_cell; ++i) {
		const auto post_connectivities = cell_elements[i].post_connectivities(post_order, connectivity_start_index);

		//vector type cast
		for (const auto& connectivity : post_connectivities) {
			const auto num_point = connectivity.size();
			std::vector<int> temp(num_point);

			for (uint j = 0; j < num_point; ++j)
				temp[j] = static_cast<int>(connectivity[j]);

			set_of_connectivities.push_back(std::move(temp));
		}

		const auto& post_nodes = set_of_post_nodes[i];
		const auto num_post_node = post_nodes.size();
		connectivity_start_index += num_post_node;
	}

	return set_of_connectivities;
}
std::vector<std::vector<double>> Grid::cell_post_coordinate_blocks(const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const {
	std::vector<std::vector<double>> coordinates(this->space_dimension_);

	for (const auto& post_nodes : set_of_post_nodes) {
		for (const auto& node : post_nodes) {
			for (ushort j = 0; j < this->space_dimension_; ++j)
				coordinates[j].push_back(node[j]);
		}
	}

	return coordinates;
}