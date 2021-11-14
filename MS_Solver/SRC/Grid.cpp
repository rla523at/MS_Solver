#include "../INC/Grid.h"

Grid::Grid(const std::string_view grid_file_path, const Grid_File_Convertor& grid_file_convertor)
{
	this->space_dimension_ = grid_file_convertor.get_space_dimension();

	auto elements = grid_file_convertor.convert_to_elements(grid_file_path);

	std::vector<Element> periodic_boundary_elements;

	for (auto& element : elements)
	{
		const auto type = element.type();
		switch (type) 
		{
		case ElementType::cell:
			this->cell_elements_.push_back(std::move(element));
			break;
		case ElementType::periodic:
			periodic_boundary_elements.push_back(std::move(element));
			break;
		default:
			this->boundary_elements_.push_back(std::move(element));
			break;
		}
	}

	this->inner_face_elements_ = this->make_inner_face_elements(periodic_boundary_elements);

}


size_t Grid::num_cells(void) const 
{
	return this->grid_elements_.cell_elements.size();
}

std::vector<Vector_Function<Polynomial>> Grid::cell_basis_vector_functions(const std::vector<ushort> solution_degrees) const 
{
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<Vector_Function<Polynomial>> basis_vector_functions(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		basis_vector_functions[i] = cell_elements[i].orthonormal_basis_vector_function(solution_degrees[i]);

	return basis_vector_functions;
}


std::vector<Euclidean_Vector> Grid::cell_center_nodes(void) const 
{
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<Euclidean_Vector> center_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		center_nodes[i] = cell_elements[i].center_node();

	return center_nodes;
}

std::vector<ushort> Grid::cell_set_of_num_post_nodes(const ushort post_order) const 
{
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

std::vector<std::vector<Euclidean_Vector>> Grid::cell_set_of_post_nodes(const ushort post_order) const 
{
	const auto& cell_elements = this->grid_elements_.cell_elements;

	const auto num_cell = cell_elements.size();
	std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		set_of_post_nodes[i] = cell_elements[i].post_nodes(post_order);

	return set_of_post_nodes;
}

std::vector<std::vector<int>> Grid::cell_set_of_connectivities(const ushort post_order, const std::vector<std::vector<Euclidean_Vector>>& set_of_post_nodes) const 
{
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

Quadrature_Rule Grid::cell_quadrature_rule(const uint cell_index, const ushort solution_degree) const
{
	const auto& cell_element = this->grid_elements_.cell_elements[cell_index];
	return cell_element.get_quadrature_rule(solution_degree);
}


//std::vector<Quadrature_Rule> Grid::cell_quadrature_rules(const std::vector<ushort> solution_degrees) const 
//{
//	const auto& cell_elements = this->grid_elements_.cell_elements;
//	const auto num_cell = cell_elements.size();
//
//	std::vector<Quadrature_Rule> quadrature_rules;
//	quadrature_rules.reserve(num_cell);
//
//	for (uint i = 0; i < num_cell; ++i)
//		quadrature_rules.push_back(cell_elements[i].get_quadrature_rule(solution_degrees[i]));
//
//	return quadrature_rules;
//}

std::vector<std::vector<double>> Grid::cell_projected_volumes(void) const
{
	const auto& cell_elements = this->grid_elements_.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<double>> projected_volumes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		projected_volumes[i] = cell_elements[i].projected_volume();

	return projected_volumes;

}

std::vector<double> Grid::cell_volumes(void) const
{
	const auto& cell_elements = this->grid_elements_.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<double> volumes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		volumes[i] = cell_elements[i].volume();

	return volumes;
}

std::vector<Element> Grid::make_inner_face_elements(const std::vector<Element>& periodic_boundary_elements)
{
	//check constructed face
	std::set<std::vector<uint>> constructed_face_vnode_index_set;

	for (const auto& boundray_element : this->boundary_elements_)
	{
		auto vnode_indexes = boundray_element.vertex_node_indexes();
		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
		constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
	}
	for (const auto& periodic_boundray_element : periodic_boundary_elements)
	{
		auto vnode_indexes = periodic_boundray_element.vertex_node_indexes();
		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
		constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
	}

	//construct face & select inner face
	std::vector<Element> inner_face_elements;

	for (const auto& cell_element : this->cell_elements_)
	{
		auto face_elements = cell_element.make_face_elements();
		for (auto& face_element : face_elements)
		{
			auto vnode_indexes = face_element.vertex_node_indexes();
			std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order

			if (constructed_face_vnode_index_set.contains(vnode_indexes))
			{
				continue;
			}
			else
			{
				constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
				inner_face_elements.push_back(std::move(face_element));
			}
		}
	}

	return inner_face_elements;
}