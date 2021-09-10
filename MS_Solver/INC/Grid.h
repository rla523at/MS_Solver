#pragma once
#include "Grid_Element_Builder.h"

#include <set>
#include <unordered_map>
#include <unordered_set>


template <ushort space_dimension>
class Grid
{
private:
	Grid_Elements<space_dimension> elements;
	std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_index_set_ignore_pbdry_;
	std::unordered_map<uint, std::set<uint>> pbdry_vnode_index_to_matched_pbdry_vnode_index_set_;
	std::unordered_map<uint, std::set<uint>> vnode_index_to_share_cell_index_set_consider_pbdry_;

public:
	Grid(Grid_Elements<space_dimension>&& grid_elements);

public:
	const Grid_Elements<space_dimension>& get_grid_elements(void) const { return this->elements; };
	const std::unordered_map<uint, std::set<uint>>& get_vnode_index_to_share_cell_index_set_consider_pbdry(void) const { return this->vnode_index_to_share_cell_index_set_consider_pbdry_; };
	const std::unordered_map<uint, std::set<uint>>& get_pbdry_vnode_index_to_matched_node_index_set(void) const { return this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_; };

public:
	std::vector<uint> boundary_owner_cell_indexes(void) const;
	std::vector<std::pair<uint, uint>> periodic_boundary_oc_nc_index_pairs(void) const;
	std::vector<std::pair<uint, uint>> inner_face_oc_nc_index_pairs(void) const;	
	std::vector<std::vector<size_t>> set_of_face_share_cell_indexes_consider_pbdry(void) const;
	std::vector<std::vector<size_t>> set_of_vertex_share_cell_indexes_consider_pbdry(void) const;

private:
	std::vector<uint> find_cell_indexes_have_these_vnodes_ignore_pbdry(const std::vector<uint>& vnode_indexes) const;
	std::vector<uint> find_cell_indexes_have_these_vnodes_consider_pbdry(const std::vector<uint>& vnode_indexes) const;
};


//Template Definition
template <ushort space_dimension>
Grid<space_dimension>::Grid(Grid_Elements<space_dimension>&& grid_elements) {
	SET_TIME_POINT;
	
	this->elements = (std::move(grid_elements));	

	//calculate vnode_index_to_share_cell_index_set_ignore_pbdry_
	const auto num_cell = this->elements.cell_elements.size();

	for (uint i = 0; i < num_cell; ++i) {
		const auto vnode_indexes = this->elements.cell_elements[i].vertex_node_indexes();
		for (const auto vnode_index : vnode_indexes) {
			if (!this->vnode_index_to_share_cell_index_set_ignore_pbdry_.contains(vnode_index))
				this->vnode_index_to_share_cell_index_set_ignore_pbdry_.emplace(vnode_index, std::set<uint>());

			this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_index).insert(i);
		}
	}

	//calculate pbdry_vnode_index_to_matched_pbdry_vnode_index_set_
	const auto num_vnode = this->vnode_index_to_share_cell_index_set_ignore_pbdry_.size();
	this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.reserve(num_vnode);

	for (const auto& [oc_side_element, nc_side_element] : this->elements.periodic_boundary_element_pairs) {
		const auto oc_side_vnode_indexes = oc_side_element.vertex_node_indexes();
		const auto nc_side_vnode_indexes = nc_side_element.vertex_node_indexes();

		const auto num_vnode = oc_side_vnode_indexes.size();
		for (ushort i = 0; i < num_vnode; ++i) {
			const auto i_vnode_index = oc_side_vnode_indexes[i];
			const auto j_vnode_index = nc_side_vnode_indexes[i];

			if (!this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.contains(i_vnode_index))
				this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.emplace(i_vnode_index, std::set<uint>());

			if (!this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.contains(j_vnode_index))
				this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.emplace(j_vnode_index, std::set<uint>());

			this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.at(i_vnode_index).insert(j_vnode_index);
			this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.at(j_vnode_index).insert(i_vnode_index);
		}


		//const auto periodic_vnode_index_pairs = oc_side_element.find_periodic_vnode_index_pairs(nc_side_element);

		//for (const auto [i_vnode_index, j_vnode_index] : periodic_vnode_index_pairs) {
			//if (!this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.contains(i_vnode_index))
			//	this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.emplace(i_vnode_index, std::set<uint>());

			//if (!this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.contains(j_vnode_index))
			//	this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.emplace(j_vnode_index, std::set<uint>());

			//this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.at(i_vnode_index).insert(j_vnode_index);
			//this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.at(j_vnode_index).insert(i_vnode_index);
		//}
	}

	//consider periodic boundary conner
	for (ushort i = 0; i < space_dimension - 1; ++i) {
		for (auto& [vnode_index, matched_vnode_index_set] : this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_) {
			if (matched_vnode_index_set.size() == 1)
				continue;

			for (const auto matched_vnode_index : matched_vnode_index_set) {
				const auto& other_matched_vnode_index_set = this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_.at(matched_vnode_index);

				auto& i_set = matched_vnode_index_set;
				const auto& j_set = other_matched_vnode_index_set;

				std::vector<uint> difference;
				std::set_difference(j_set.begin(), j_set.end(), i_set.begin(), i_set.end(), std::back_inserter(difference));

				i_set.insert(difference.begin(), difference.end());
				i_set.erase(vnode_index);
			}
		}
	}

	//calculate vnode_index_to_share_cell_index_set_consider_pbdry_
	this->vnode_index_to_share_cell_index_set_consider_pbdry_= this->vnode_index_to_share_cell_index_set_ignore_pbdry_;

	for (const auto& [vnode_index, matched_vnode_index_set] : this->pbdry_vnode_index_to_matched_pbdry_vnode_index_set_) {
		for (const auto matched_vnode_index : matched_vnode_index_set) {
			auto& i_set = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_index);
			auto& j_set = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(matched_vnode_index);

			std::vector<uint> difference;
			std::set_difference(i_set.begin(), i_set.end(), j_set.begin(), j_set.end(), std::back_inserter(difference));

			i_set.insert(difference.begin(), difference.end());
			j_set.insert(difference.begin(), difference.end());
		}
	}

	Log::content_ << std::left << std::setw(50) << "@ Make grid " << " ----------- " << GET_TIME_DURATION << "s\n\n";
	Log::print();
}

template <ushort space_dimension>
std::vector<uint> Grid<space_dimension>::boundary_owner_cell_indexes(void) const {			
	const auto& boundary_elements = this->elements.boundary_elements;

	const auto num_boundary = boundary_elements.size();
	std::vector<uint> boundary_oc_indexes(num_boundary);

	for (uint i = 0; i < num_boundary; ++i) {
		const auto& boundary_element = boundary_elements[i];

		const auto vnode_indexes = boundary_element.vertex_node_indexes();
		const auto oc_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(vnode_indexes);
		dynamic_require(oc_indexes.size() == 1, "boundary should have unique owner cell");

		const auto oc_index = oc_indexes.front();
		boundary_oc_indexes[i] = oc_index;
	}

	return boundary_oc_indexes;
}

template <ushort space_dimension>
std::vector<std::pair<uint, uint>> Grid<space_dimension>::periodic_boundary_oc_nc_index_pairs(void) const {
	const auto& pbdry_element_pairs = this->elements.periodic_boundary_element_pairs;

	const auto num_pbdry_pair = pbdry_element_pairs.size();
	std::vector<std::pair<uint, uint>> pbdry_oc_nc_index_pairs(num_pbdry_pair);

	for (size_t i = 0; i < num_pbdry_pair; ++i) {
		const auto& [oc_side_element, nc_side_element] = pbdry_element_pairs[i];

		const auto oc_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(oc_side_element.vertex_node_indexes());
		const auto nc_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(nc_side_element.vertex_node_indexes());
		dynamic_require(oc_indexes.size() == 1, "periodic boundary should have unique owner cell");
		dynamic_require(nc_indexes.size() == 1, "periodic boundary should have unique neighbor cell");

		const auto oc_index = oc_indexes.front();
		const auto nc_index = nc_indexes.front();
		pbdry_oc_nc_index_pairs[i] = { oc_index,nc_index };
	}

	return pbdry_oc_nc_index_pairs;
}

template <ushort space_dimension>
std::vector<std::pair<uint, uint>> Grid<space_dimension>::inner_face_oc_nc_index_pairs(void) const {
	const auto& inner_face_elements = this->elements.inner_face_elements;

	const auto num_inner_face = inner_face_elements.size();
	std::vector<std::pair<uint, uint>> inner_face_oc_nc_index_pairs(num_inner_face);

	for (uint i = 0; i < num_inner_face; ++i) {
		const auto& inner_face_element = inner_face_elements[i];

		const auto cell_indexes = this->find_cell_indexes_have_these_vnodes_ignore_pbdry(inner_face_element.vertex_node_indexes());
		dynamic_require(cell_indexes.size() == 2, "inner face should have an unique owner neighbor cell pair");

		//set first index as oc index
		const auto oc_index = cell_indexes[0];
		const auto nc_index = cell_indexes[1];
		inner_face_oc_nc_index_pairs[i] = { oc_index,nc_index };
	}

	return inner_face_oc_nc_index_pairs;
}

template <ushort space_dimension>
std::vector<std::vector<size_t>> Grid<space_dimension>::set_of_face_share_cell_indexes_consider_pbdry(void) const {
	const auto& cell_elements = this->elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<size_t>> set_of_face_share_cell_indexes;
	set_of_face_share_cell_indexes.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;

		const auto face_vnode_indexes_set = element.set_of_face_vertex_node_indexes();
		const auto num_face = face_vnode_indexes_set.size();

		std::vector<size_t> face_share_cell_indexes;
		face_share_cell_indexes.reserve(num_face);

		for (const auto& face_vnode_indexes : face_vnode_indexes_set) {
			auto this_face_share_cell_indexes = this->find_cell_indexes_have_these_vnodes_consider_pbdry(face_vnode_indexes);

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
std::vector<std::vector<size_t>> Grid<space_dimension>::set_of_vertex_share_cell_indexes_consider_pbdry(void) const {
	const auto& cell_elements = this->elements.cell_elements;
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<size_t>> set_of_vertex_share_cell_indexes;
	set_of_vertex_share_cell_indexes.reserve(num_cell);

	for (size_t i = 0; i < num_cell; ++i) {
		const auto& element = cell_elements[i];
		const auto& geometry = element.geometry_;
				
		auto vnode_indexes = element.vertex_node_indexes();
		std::set<size_t> vertex_share_cell_index_set;
		for (const auto vnode_index : vnode_indexes) {
			const auto& vnode_share_cell_index_set = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_index);
			vertex_share_cell_index_set.insert(vnode_share_cell_index_set.begin(), vnode_share_cell_index_set.end());
		}
		vertex_share_cell_index_set.erase(i);

		set_of_vertex_share_cell_indexes.push_back({ vertex_share_cell_index_set.begin(), vertex_share_cell_index_set.end() });
	}

	return set_of_vertex_share_cell_indexes;
}




template<ushort space_dimension>
std::vector<uint> Grid<space_dimension>::find_cell_indexes_have_these_vnodes_ignore_pbdry(const std::vector<uint>& vnode_indexes) const {
	std::vector<uint> this_face_share_cell_indexes;

	const auto& set0 = this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_indexes[0]);
	const auto& set1 = this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_indexes[1]);
	std::set_intersection(set0.begin(), set0.end(), set1.begin(), set1.end(), std::back_inserter(this_face_share_cell_indexes));

	const auto num_face_vnode = vnode_indexes.size();

	if (2 < num_face_vnode) {
		std::vector<uint> temp;
		for (ushort i = 2; i < num_face_vnode; ++i) {
			const auto& set_i = this->vnode_index_to_share_cell_index_set_ignore_pbdry_.at(vnode_indexes[i]);

			std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(temp));
			std::swap(this_face_share_cell_indexes, temp);
			temp.clear();
		}
	}

	return this_face_share_cell_indexes;
}

template<ushort space_dimension>
std::vector<uint> Grid<space_dimension>::find_cell_indexes_have_these_vnodes_consider_pbdry(const std::vector<uint>& vnode_indexes) const {
	std::vector<uint> this_face_share_cell_indexes;

	const auto& set0 = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_indexes[0]);
	const auto& set1 = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_indexes[1]);
	std::set_intersection(set0.begin(), set0.end(), set1.begin(), set1.end(), std::back_inserter(this_face_share_cell_indexes));

	const auto num_face_vnode = vnode_indexes.size();

	if (2 < num_face_vnode) {
		std::vector<uint> temp;
		for (ushort i = 2; i < num_face_vnode; ++i) {
			const auto& set_i = this->vnode_index_to_share_cell_index_set_consider_pbdry_.at(vnode_indexes[i]);

			std::set_intersection(this_face_share_cell_indexes.begin(), this_face_share_cell_indexes.end(), set_i.begin(), set_i.end(), std::back_inserter(temp));
			std::swap(this_face_share_cell_indexes, temp);
			temp.clear();
		}
	}

	return this_face_share_cell_indexes;
}