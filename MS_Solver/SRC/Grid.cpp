#include "../INC/Grid.h"

Grid::Grid(const std::string_view grid_file_path, const Grid_File_Convertor& grid_file_convertor)
{
	this->space_dimension_ = grid_file_convertor.get_space_dimension();
	auto elements = grid_file_convertor.convert_to_elements(grid_file_path);

	Profiler::set_time_point();

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
	
	this->periodic_boundary_element_pairs_ = this->make_periodic_boundary_element_pairs(std::move(periodic_boundary_elements));
	this->inner_face_elements_ = this->make_inner_face_elements();

	LOG << std::left << std::setw(50) << "@ Make Grid Element" << " ----------- " << Profiler::get_time_duration() << "s\n\n";
	LOG << "  " << std::setw(8) << this->cell_elements_.size() << " cell \n";
	LOG << "  " << std::setw(8) << this->boundary_elements_.size() << " boundary\n";
	LOG << "  " << std::setw(8) << this->periodic_boundary_element_pairs_.size() << " periodic boundary pair\n";
	LOG << "  " << std::setw(8) << this->inner_face_elements_.size() << " inner face\n\n" << Log::print_;


	Profiler::set_time_point();

	//calculate vnode_index_to_share_cell_index_set_ignore_pbdry_
	const auto num_cell = this->cell_elements_.size();

	for (uint i = 0; i < num_cell; ++i) 
	{
		const auto vnode_indexes = this->cell_elements_[i].vertex_node_indexes();
		for (const auto vnode_index : vnode_indexes) 
		{
			if (!this->vnode_index_to_share_cell_index_set_ignore_pbdry_.contains(vnode_index))
			{
				this->vnode_index_to_share_cell_index_set_ignore_pbdry_.emplace(vnode_index, std::set<uint>());
			}

			this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_index).insert(i);
		}
	}

	//calculate vnode_index_to_share_cell_index_set_consider_pbdry_
	this->vnode_index_to_share_cell_index_set_consider_pbdry_ = this->vnode_index_to_share_cell_index_set_ignore_pbdry_;
	
	const auto vnode_index_to_periodic_matched_node_index_set = this->calculate_vnode_index_to_peridoic_matched_node_index_set();

	for (const auto& [vnode_index, periodic_matched_node_index_set] : vnode_index_to_periodic_matched_node_index_set)
	{
		for (const auto matched_vnode_index : periodic_matched_node_index_set) 
		{
			auto& i_set = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_index);
			const auto& j_set = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(matched_vnode_index);

			const auto difference = ms::set_difference(j_set, i_set);

			i_set.insert(difference.begin(), difference.end());
		}
	}

	LOG << std::left << std::setw(50) << "@ Make grid connecitivy " << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;	
}

ushort Grid::space_dimension(void) const
{
	return this->space_dimension_;
}


size_t Grid::num_cells(void) const 
{
	return this->cell_elements_.size();
}

std::vector<Vector_Function<Polynomial>> Grid::cell_basis_vector_functions(const std::vector<ushort> solution_degrees) const 
{
	const auto& cell_elements = this->cell_elements_;

	const auto num_cell = cell_elements.size();
	std::vector<Vector_Function<Polynomial>> basis_vector_functions(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		basis_vector_functions[i] = cell_elements[i].orthonormal_basis_vector_function(solution_degrees[i]);

	return basis_vector_functions;
}


std::vector<Euclidean_Vector> Grid::cell_center_nodes(void) const 
{
	const auto& cell_elements = this->cell_elements_;

	const auto num_cell = cell_elements.size();
	std::vector<Euclidean_Vector> center_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		center_nodes[i] = cell_elements[i].center_node();

	return center_nodes;
}

std::vector<ushort> Grid::cell_set_of_num_post_points(const ushort post_order) const 
{
	const auto& cell_elements = this->cell_elements_;

	const auto num_cell = cell_elements.size();
	std::vector<ushort> cell_set_of_num_post_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
	{
		cell_set_of_num_post_nodes[i] = cell_elements[i].num_post_nodes(post_order);
	}

	return cell_set_of_num_post_nodes;
}

std::vector<ushort> Grid::cell_set_of_num_post_elements(const ushort post_order) const{
	const auto& cell_elements = this->cell_elements_;

	const auto num_cell = cell_elements.size();
	std::vector<ushort> cell_set_of_num_post_elements(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		cell_set_of_num_post_elements[i] = cell_elements[i].num_post_elements(post_order);

	return cell_set_of_num_post_elements;
}

std::vector<std::vector<Euclidean_Vector>> Grid::cell_set_of_post_element_centers(const ushort post_order) const
{
	const auto num_cell = this->cell_elements_.size();
	std::vector<std::vector<Euclidean_Vector>> set_of_post_element_centers(num_cell);

	for (int i = 0; i < num_cell; ++i)
	{
		set_of_post_element_centers[i] = this->cell_elements_[i].post_element_centers(post_order);
	}

	return set_of_post_element_centers;
}

std::vector<std::vector<Euclidean_Vector>> Grid::cell_set_of_post_points(const ushort post_order) const 
{
	const auto num_cell = this->cell_elements_.size();
	std::vector<std::vector<Euclidean_Vector>> set_of_post_nodes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
	{
		set_of_post_nodes[i] = this->cell_elements_[i].post_nodes(post_order);
	}

	return set_of_post_nodes;
}

std::vector<std::vector<int>> Grid::cell_set_of_connectivities(const ushort post_order, const std::vector<std::vector<Euclidean_Vector>>& set_of_post_points) const 
{
	const auto& cell_elements = this->cell_elements_;

	size_t connectivity_start_index = 0;

	std::vector<std::vector<int>> set_of_connectivities;

	for (uint i = 0; i < cell_elements.size(); ++i) 
	{
		auto icell_connectivities = cell_elements[i].post_connectivities(post_order, connectivity_start_index);
		ms::merge(set_of_connectivities, std::move(icell_connectivities));

		const auto& post_points = set_of_post_points[i];
		const auto num_post_points = post_points.size();
		connectivity_start_index += num_post_points;
	}

	return set_of_connectivities;
}

Quadrature_Rule Grid::cell_quadrature_rule(const uint cell_index, const ushort solution_degree) const
{
	const auto& cell_element = this->cell_elements_[cell_index];
	return cell_element.get_quadrature_rule(solution_degree);
}

std::vector<std::vector<double>> Grid::cell_projected_volumes(void) const
{
	const auto& cell_elements = this->cell_elements_;
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<double>> projected_volumes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		projected_volumes[i] = cell_elements[i].projected_volume();

	return projected_volumes;

}

std::vector<double> Grid::cell_volumes(void) const
{
	const auto& cell_elements = this->cell_elements_;
	const auto num_cell = cell_elements.size();

	std::vector<double> volumes(num_cell);

	for (uint i = 0; i < num_cell; ++i)
		volumes[i] = cell_elements[i].volume();

	return volumes;
}


size_t Grid::num_boundaries(void) const
{
	return this->boundary_elements_.size();
}

uint Grid::boundary_owner_cell_index(const uint bdry_index) const
{
	const auto vnode_indexes = this->boundary_elements_[bdry_index].vertex_node_indexes();
	const auto cell_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(vnode_indexes);
	REQUIRE(cell_indexes.size() == 1, "boundary should have unique owner cell");

	return cell_indexes.front();
}

std::vector<uint> Grid::boundary_owner_cell_indexes(void) const
{
	const auto num_boundary = this->boundary_elements_.size();
	std::vector<uint> oc_indexes(num_boundary);

	for (uint i = 0; i < num_boundary; ++i)
	{
		oc_indexes[i] = this->boundary_owner_cell_index(i);
	}

	return oc_indexes;
}

ElementType Grid::boundary_type(const uint bdry_index) const
{
	return this->boundary_elements_[bdry_index].type();
}

std::vector<ElementType> Grid::boundary_types(void) const {
	const auto num_boundary = this->boundary_elements_.size();
	std::vector<ElementType> types(num_boundary);

	for (uint i = 0; i < num_boundary; ++i)
	{
		types[i] = this->boundary_type(i);
	}

	return types;
}

Quadrature_Rule Grid::boundary_quadrature_rule(const uint bdry_index, const ushort polynomial_degree) const
{
	return this->boundary_elements_[bdry_index].get_quadrature_rule(polynomial_degree);
}

std::vector<Euclidean_Vector> Grid::boundary_normals(const uint bdry_index, const uint oc_indexes, const std::vector<Euclidean_Vector>& points) const
{
	const auto& bdry_element = this->boundary_elements_[bdry_index];
	const auto& oc_element = this->cell_elements_[oc_indexes];

	return bdry_element.outward_normalized_normal_vectors(oc_element, points);
}

size_t Grid::num_inner_faces(void) const
{
	return this->inner_face_elements_.size();
}

std::pair<uint, uint> Grid::inner_face_oc_nc_index_pair(const uint inner_face_index) const
{
	const auto cell_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(this->inner_face_elements_[inner_face_index].vertex_node_indexes());
	REQUIRE(cell_indexes.size() == 2, "inner face should have an unique owner neighbor cell pair");

	//set first index as oc index
	const auto oc_index = cell_indexes[0];
	const auto nc_index = cell_indexes[1];
	return { oc_index,nc_index };
}

Quadrature_Rule Grid::inner_face_quadrature_rule(const uint inner_face_index, const ushort polynomial_degree) const
{
	return  this->inner_face_elements_[inner_face_index].get_quadrature_rule(polynomial_degree);
}

std::vector<Euclidean_Vector> Grid::inner_face_normals(const uint inner_face_index, const uint oc_index, const std::vector<Euclidean_Vector>& points) const
{
	const auto& inner_face_element = this->inner_face_elements_[inner_face_index];
	const auto& oc_element = this->cell_elements_[oc_index];

	return inner_face_element.outward_normalized_normal_vectors(oc_element, points);
}

size_t Grid::num_periodic_boundary_pairs(void) const
{
	return this->periodic_boundary_element_pairs_.size();
}

std::pair<uint, uint> Grid::periodic_boundary_oc_nc_index_pair(const uint pbdry_pair_index) const
{
	const auto& [oc_side_element, nc_side_element] = this->periodic_boundary_element_pairs_[pbdry_pair_index];

	const auto oc_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(oc_side_element.vertex_node_indexes());
	const auto nc_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(nc_side_element.vertex_node_indexes());
	REQUIRE(oc_indexes.size() == 1, "periodic boundary should have unique owner cell");
	REQUIRE(nc_indexes.size() == 1, "periodic boundary should have unique neighbor cell");

	const auto oc_index = oc_indexes.front();
	const auto nc_index = nc_indexes.front();
	return { oc_index,nc_index };
}

std::pair<Quadrature_Rule, Quadrature_Rule> Grid::periodic_boundary_quadrature_rule_pair(const uint pbdry_pair_index, const ushort solution_degree) const
{
	const auto& [oc_side_element, nc_side_element] = this->periodic_boundary_element_pairs_[pbdry_pair_index];
	return { oc_side_element.get_quadrature_rule(solution_degree), nc_side_element.get_quadrature_rule(solution_degree) };
}

std::vector<uint> Grid::find_cell_indexes_have_these_vnodes_ignore_pbdry(const std::vector<uint>& vnode_indexes) const
{
	const auto num_vnode = vnode_indexes.size();

	std::vector<const std::set<uint>*> share_cell_index_set_ptrs;
	share_cell_index_set_ptrs.reserve(num_vnode);

	for (int i = 0; i < num_vnode; ++i)
	{
		share_cell_index_set_ptrs.push_back(&this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_indexes[i]));
	}

	return ms::set_intersection(share_cell_index_set_ptrs);
}

std::vector<Element> Grid::make_inner_face_elements(void) const
{
	std::set<std::vector<uint>> constructed_face_vnode_index_set;

	for (const auto& boundray_element : this->boundary_elements_)
	{
		auto vnode_indexes = boundray_element.vertex_node_indexes();
		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
		constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
	}

	for (const auto& [pbdry_elem1, pbdry_elem2] : this->periodic_boundary_element_pairs_)
	{
		auto vnode_indexes1 = pbdry_elem1.vertex_node_indexes();
		auto vnode_indexes2 = pbdry_elem2.vertex_node_indexes();

		std::sort(vnode_indexes1.begin(), vnode_indexes1.end());	//to ignore index order
		std::sort(vnode_indexes2.begin(), vnode_indexes2.end());	//to ignore index order

		constructed_face_vnode_index_set.insert(std::move(vnode_indexes1));
		constructed_face_vnode_index_set.insert(std::move(vnode_indexes2));
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

std::vector<std::pair<Element, Element>> Grid::make_periodic_boundary_element_pairs(std::vector<Element>&& periodic_boundary_elements) const
{
	const auto num_periodic_element = periodic_boundary_elements.size();
	const auto num_pair = static_cast<uint>(0.5 * num_periodic_element);

	std::unordered_set<uint> matched_index_set;
	matched_index_set.reserve(num_periodic_element);

	std::vector<std::pair<Element, Element>> matched_periodic_element_pairs;
	matched_periodic_element_pairs.reserve(num_pair);

	for (uint i = 0; i < num_periodic_element; ++i)
	{
		if (matched_index_set.contains(i))
		{
			continue;
		}

		auto& i_element = periodic_boundary_elements[i];

		for (uint j = i + 1; j < num_periodic_element; ++j)
		{
			if (matched_index_set.contains(j))
			{
				continue;
			}

			auto& j_element = periodic_boundary_elements[j];

			auto periodic_matched_node_indexes = j_element.find_periodic_matched_node_indexes(i_element);

			if (!periodic_matched_node_indexes.empty())
			{
				j_element.rearrange_node_indexes(std::move(periodic_matched_node_indexes));

				matched_periodic_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(j_element)));
				matched_index_set.insert(i);
				matched_index_set.insert(j);
				break;
			}
		}
	}

	REQUIRE(matched_periodic_element_pairs.size() == num_pair, "every periodic boundary should be matched");
	return matched_periodic_element_pairs;
}

std::unordered_map<uint, std::set<uint>> Grid::calculate_vnode_index_to_peridoic_matched_node_index_set(void) const
{
	std::unordered_map<uint, std::set<uint>> vnode_index_to_periodic_matched_vnode_index_set;

	for (const auto& [oc_side_element, nc_side_element] : this->periodic_boundary_element_pairs_) 
	{
		const auto oc_side_vnode_indexes = oc_side_element.vertex_node_indexes();
		const auto nc_side_vnode_indexes = nc_side_element.vertex_node_indexes();

		const auto num_vnode = oc_side_vnode_indexes.size();
		for (ushort i = 0; i < num_vnode; ++i) 
		{
			const auto i_vnode_index = oc_side_vnode_indexes[i];
			const auto j_vnode_index = nc_side_vnode_indexes[i];

			if (!vnode_index_to_periodic_matched_vnode_index_set.contains(i_vnode_index))
				vnode_index_to_periodic_matched_vnode_index_set.emplace(i_vnode_index, std::set<uint>());

			if (!vnode_index_to_periodic_matched_vnode_index_set.contains(j_vnode_index))
				vnode_index_to_periodic_matched_vnode_index_set.emplace(j_vnode_index, std::set<uint>());

			vnode_index_to_periodic_matched_vnode_index_set.at(i_vnode_index).insert(j_vnode_index);
			vnode_index_to_periodic_matched_vnode_index_set.at(j_vnode_index).insert(i_vnode_index);
		}
	}

	//consider pbdry conner
	for (ushort i = 0; i < this->space_dimension_ - 1; ++i) 
	{
		for (auto& [vnode_index, matched_vnode_index_set] : vnode_index_to_periodic_matched_vnode_index_set) 
		{
			if (matched_vnode_index_set.size() == 1)
			{
				continue;
			}

			for (const auto matched_vnode_index : matched_vnode_index_set) 
			{
				const auto& other_matched_vnode_index_set = vnode_index_to_periodic_matched_vnode_index_set.at(matched_vnode_index);

				auto& i_set = matched_vnode_index_set;
				const auto& j_set = other_matched_vnode_index_set;

				const auto difference = ms::set_difference(j_set, i_set);

				if (difference.empty())
				{
					continue;
				}

				i_set.insert(difference.begin(), difference.end());
				i_set.erase(vnode_index);
			}
		}
	}

	return vnode_index_to_periodic_matched_vnode_index_set;
}