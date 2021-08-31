#pragma once
#include "Grid_Element_Builder.h"

#include <set>
#include <unordered_map>
#include <unordered_set>


struct Grid_Connectivity
{
	std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_index_set; // vnode				:= node at vertex, 
																				// share_cell_indexes	:= cell indexes which shares vnode
	std::vector<uint> boundary_oc_indexes;

	std::vector<std::pair<uint, uint>> periodic_boundary_oc_nc_index_pairs;			// {owner cell index, neighbor cell index}
	std::unordered_map<uint, std::set<uint>> vnode_index_to_matched_vnode_index_set;

	std::vector<std::pair<uint, uint>> inner_face_oc_nc_index_pairs;				// {owner cell index, neighbor cell index}
};


template <ushort space_dimension>
struct Grid
{
	Grid_Elements<space_dimension> elements;
	Grid_Connectivity connectivity;
		
	std::vector<std::vector<size_t>> calculate_set_of_face_share_cell_indexes(void) const;
	std::vector<std::vector<size_t>> calculate_set_of_vertex_share_cell_indexes(void) const;
};


template<ushort space_dimension>
class Grid_Builder
{
private:
	Grid_Builder(void) = delete;

private:
	using This_				= Grid_Builder;
	using Space_Vector_		= Euclidean_Vector<space_dimension>;

public:
	template <typename Grid_File_Type>
	static Grid<space_dimension> build(const std::string& grid_file_name);

private:
	static Grid_Connectivity make_grid_connectivity(const Grid_Elements<space_dimension>& grid_elements);
	static std::vector<uint> find_cell_indexes_have_these_vnodes(const std::unordered_map<uint, std::set<uint>>& vertex_node_index_to_cell_container_indexes, const std::vector<uint>& face_node_indexes);
};


//template definition part
template <ushort space_dimension>
std::vector<std::vector<size_t>> Grid<space_dimension>::calculate_set_of_face_share_cell_indexes(void) const {
	const auto& vnode_index_to_share_cell_index_set = this->connectivity.vnode_index_to_share_cell_index_set;
	const auto& cell_elements = this->elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<size_t>> set_of_face_share_cell_indexes;
	set_of_face_share_cell_indexes.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;

		const auto face_vnode_indexes_set = element.face_vertex_node_indexes_set();
		const auto num_face = face_vnode_indexes_set.size();

		std::vector<size_t> face_share_cell_indexes;
		face_share_cell_indexes.reserve(num_face);

		for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
			std::vector<size_t> this_face_share_cell_indexes;

			const auto num_face_vnode = face_vnode_indexes.size();

			const auto& set_0 = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[0]);
			const auto& set_1 = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[1]);
			std::set_intersection(set_0.begin(), set_0.end(), set_1.begin(), set_1.end(), std::back_inserter(this_face_share_cell_indexes));

			if (2 < num_face_vnode) {
				std::vector<size_t> temp;
				for (size_t i = 2; i < num_face_vnode; ++i) {
					const auto& set_i = vnode_index_to_share_cell_index_set.at(face_vnode_indexes[i]);

					std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(temp));
					std::swap(this_face_share_cell_indexes, temp);
					temp.clear();
				}
			}

			const auto my_index_pos_iter = std::find(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), i);
			dynamic_require(my_index_pos_iter != this_face_share_cell_indexes.end(), "my index should be included in this face share cell indexes");

			this_face_share_cell_indexes.erase(my_index_pos_iter);
			dynamic_require(this_face_share_cell_indexes.size() <= 1, "face share cell should be unique or absent");

			face_share_cell_indexes.push_back(this_face_share_cell_indexes.front());
		}

		set_of_face_share_cell_indexes.push_back(std::move(face_share_cell_indexes));
	}

	return set_of_face_share_cell_indexes; //this is not ordered set
}


template <ushort space_dimension>
std::vector<std::vector<size_t>> Grid<space_dimension>::calculate_set_of_vertex_share_cell_indexes(void) const {
	const auto& cell_elements = this->elements.cell_elements;
	const auto& vnode_index_to_share_cell_index_set = this->connectivity.vnode_index_to_share_cell_index_set;

	const auto num_cell = cell_elements.size();
	std::vector<std::vector<size_t>> set_of_vertex_share_cell_indexes;
	set_of_vertex_share_cell_indexes.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;
				
		auto vnode_indexes = element.vertex_node_indexes();
		std::set<size_t> vertex_share_cell_index_set;
		for (const auto vnode_index : vnode_indexes) {
			const auto& share_cell_indexes = vnode_index_to_share_cell_index_set.at(vnode_index);
			vertex_share_cell_index_set.insert(share_cell_indexes.begin(), share_cell_indexes.end());
		}
		vertex_share_cell_index_set.erase(i);

		set_of_vertex_share_cell_indexes.push_back({ vertex_share_cell_index_set.begin(), vertex_share_cell_index_set.end() });
	}

	return set_of_vertex_share_cell_indexes;
}

template <ushort space_dimension>
template <typename Grid_File_Type>
Grid<space_dimension> Grid_Builder<space_dimension>::build(const std::string& grid_file_name) {
	static_require(ms::is_grid_file_type<Grid_File_Type>, "It should be grid file type");

	const auto grid_elements		= Grid_Element_Builder<Grid_File_Type, space_dimension>::build_from_grid_file(grid_file_name);
	const auto grid_connectivity	= This_::make_grid_connectivity(grid_elements);
	return { grid_elements,grid_connectivity };
}


template <ushort space_dimension>
Grid_Connectivity Grid_Builder<space_dimension>::make_grid_connectivity(const Grid_Elements<space_dimension>& grid_elements) {
	SET_TIME_POINT;

	const auto& [cell_elements, boundary_elements, pbdry_element_pairs, inner_face_elements] = grid_elements;

	std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_index_set;

	const auto num_cell = cell_elements.size();
	for (uint i = 0; i < num_cell; ++i) {
		const auto vnode_indexes = cell_elements[i].vertex_node_indexes();
		for (const auto& vnode_index : vnode_indexes) {
			if (!vnode_index_to_share_cell_index_set.contains(vnode_index))
				vnode_index_to_share_cell_index_set.emplace(vnode_index, std::set<uint>());

			vnode_index_to_share_cell_index_set.at(vnode_index).insert(i);
		}
	}

	//boudnary grid connectivity
	const auto num_boundary = boundary_elements.size();
	std::vector<uint> boudnary_oc_indexes(num_boundary);

	for (uint i = 0; i < num_boundary; ++i) {
		const auto& boundary_element = boundary_elements[i];

		const auto vnode_indexes = boundary_element.vertex_node_indexes();
		const auto oc_indexes = This_::find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_index_set, vnode_indexes);
		dynamic_require(oc_indexes.size() == 1, "boundary should have unique owner cell");

		const auto oc_index = oc_indexes.front();
		boudnary_oc_indexes[i] = oc_index;
	}

	//periodic boundary grid connectivity
	const auto num_vnode = vnode_index_to_share_cell_index_set.size();
	const auto num_pbdry_pair = pbdry_element_pairs.size();
	std::vector<std::pair<uint, uint>> pbdry_oc_nc_index_pairs(num_pbdry_pair);	
	std::unordered_map<uint, std::set<uint>> vnode_index_to_matched_vnode_index_set;	
	vnode_index_to_matched_vnode_index_set.reserve(num_vnode);

	for (size_t i = 0; i < num_pbdry_pair; ++i) {
		const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];

		const auto periodic_vnode_index_pairs = oc_side_element.find_periodic_vnode_index_pairs(nc_side_element);
		for (const auto [i_vnode_index, j_vnode_index] : periodic_vnode_index_pairs) {
			if (!vnode_index_to_matched_vnode_index_set.contains(i_vnode_index))
				vnode_index_to_matched_vnode_index_set.emplace(i_vnode_index, std::set<uint>());

			if (!vnode_index_to_matched_vnode_index_set.contains(j_vnode_index))
				vnode_index_to_matched_vnode_index_set.emplace(j_vnode_index, std::set<uint>());

			vnode_index_to_matched_vnode_index_set.at(i_vnode_index).insert(j_vnode_index);
			vnode_index_to_matched_vnode_index_set.at(j_vnode_index).insert(i_vnode_index);
		}

		const auto oc_indexes = This_::find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_index_set, oc_side_element.vertex_node_indexes());
		const auto nc_indexes = This_::find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_index_set, nc_side_element.vertex_node_indexes());
		dynamic_require(oc_indexes.size() == 1, "periodic boundary should have unique owner cell");
		dynamic_require(nc_indexes.size() == 1, "periodic boundary should have unique neighbor cell");

		const auto oc_index = oc_indexes.front();
		const auto nc_index = nc_indexes.front();
		pbdry_oc_nc_index_pairs[i] = { oc_index,nc_index };
	}

	// check missing match - periodic boundary conner
	for (ushort i = 0; i < space_dimension - 1; ++i) {
		for (auto& [vnode_index, matched_vnode_index_set] : vnode_index_to_matched_vnode_index_set) {
			if (matched_vnode_index_set.size() == 1)
				continue;

			for (const auto matched_vnode_index : matched_vnode_index_set) {
				const auto& other_matched_vnode_index_set = vnode_index_to_matched_vnode_index_set.at(matched_vnode_index);

				//other matched vnode index should be included my matched vnode index
				for (const auto other_matched_vnode_index : other_matched_vnode_index_set) {
					if (other_matched_vnode_index != vnode_index) //except my vnode index
						matched_vnode_index_set.insert(other_matched_vnode_index);
				}
			}
		}
	}

	// reflect pbdry connectivity to vnode_index_to_share_cell_indexes
	for (const auto& [vnode_index, matched_vnode_index_set] : vnode_index_to_matched_vnode_index_set) {
		for (const auto matched_vnode_index : matched_vnode_index_set) {
			auto& i_set = vnode_index_to_share_cell_index_set.at(vnode_index);
			auto& j_set = vnode_index_to_share_cell_index_set.at(matched_vnode_index);

			//union
			i_set.insert(j_set.begin(), j_set.end());
			j_set.insert(i_set.begin(), i_set.end());
		}
	}

	//inner face grid connectivity
	const auto num_inner_face = inner_face_elements.size();
	std::vector<std::pair<uint, uint>> inner_face_oc_nc_index_pairs(num_inner_face);

	for (uint i = 0; i < num_inner_face; ++i) {
		const auto& inner_face_element = inner_face_elements[i];

		const auto cell_indexes = This_::find_cell_indexes_have_these_vnodes(vnode_index_to_share_cell_index_set, inner_face_element.vertex_node_indexes());
		dynamic_require(cell_indexes.size() == 2, "inner face should have an unique owner neighbor cell pair");

		//set first index as oc index
		const auto oc_index = cell_indexes[0];
		const auto nc_index = cell_indexes[1];
		inner_face_oc_nc_index_pairs[i] = { oc_index,nc_index };
	}

	Log::content_ << std::left << std::setw(50) << "@ Figure out connectivity" << " ----------- " << GET_TIME_DURATION << "s\n\n";
	Log::print();

	return { vnode_index_to_share_cell_index_set,	boudnary_oc_indexes, pbdry_oc_nc_index_pairs, vnode_index_to_matched_vnode_index_set, inner_face_oc_nc_index_pairs };
}


template<ushort space_dimension>
std::vector<uint> Grid_Builder<space_dimension>::find_cell_indexes_have_these_vnodes(const std::unordered_map<uint, std::set<uint>>& vertex_node_index_to_cell_container_indexes, const std::vector<uint>& face_node_indexes) {
	const auto start_node_index = face_node_indexes[0];
	const auto end_node_index = face_node_indexes[1];

	const auto& indexes_have_start_node = vertex_node_index_to_cell_container_indexes.at(start_node_index);
	const auto& indexes_have_end_node = vertex_node_index_to_cell_container_indexes.at(end_node_index);

	std::vector<uint> cell_continaer_indexes_have_these_nodes;
	std::set_intersection(indexes_have_start_node.begin(), indexes_have_start_node.end(), indexes_have_end_node.begin(), indexes_have_end_node.end(), std::back_inserter(cell_continaer_indexes_have_these_nodes));

	return cell_continaer_indexes_have_these_nodes;
}