#pragma once
#include "Geometry.h"

#include <unordered_set>

enum class ElementType
{
	cell, face,
	slip_wall,
	supersonic_inlet1, supersonic_inlet2,
	supersonic_outlet,
	initial_constant_BC,
	periodic,
	not_in_list
};

enum class FaceType
{
	inward_face,
	outward_face,
	not_my_face
};

class Element : public Geometry
{
public:
	Element(const ElementType element_type, std::vector<uint>&& node_indexes, Geometry&& geometry)
		: element_type_(element_type), 
		node_indexes_(std::move(node_indexes)),
		Geometry(std::move(geometry)) {};

public://Query
	std::vector<uint> find_periodic_match_node_sequences(const Element& other) const;
	bool is_periodic_pair(const Element& other) const;
	std::vector<Element> make_face_elements(void) const;
	Euclidean_Vector outward_normalized_normal_vector(const Element& owner_cell_element, const Euclidean_Vector& point) const;
	ElementType type(void) const;
	std::vector<uint> vertex_node_indexes(void) const;

private:
	FaceType check_face_type(const Element& owner_cell_element) const;
	bool is_periodic_boundary(void) const;
	std::vector<std::vector<uint>> set_of_face_node_indexes(void) const;
	std::vector<std::vector<uint>> set_of_face_vertex_node_indexes(void) const;


private:
	ElementType element_type_;
	std::vector<uint> node_indexes_;


	//std::vector<Euclidean_Vector> normalized_normal_vectors(const Element& owner_cell_element, const std::vector<Euclidean_Vector>& nodes) const;
	
	//Euclidean_Vector node_at_index(const uint index) const;
	//std::vector<Euclidean_Vector> nodes_at_indexes(const std::vector<uint>& indexes) const;


	////private:

};
















//std::vector<Euclidean_Vector> Element::normalized_normal_vectors(const Element& owner_cell_element, const std::vector<Euclidean_Vector>& nodes) const {
//	const auto num_node = nodes.size();	
//	std::vector<Euclidean_Vector> normalized_normal_vectors(num_node);
//
//	for (uint i = 0; i < num_node; ++i)
//		normalized_normal_vectors[i] = this->geometry_.normalized_normal_vector(nodes[i]);
//
//	const auto face_type = this->check_face_type(owner_cell_element);
//	if (face_type == FaceType::inward_face) {
//		for (auto& normal_vector : normalized_normal_vectors)
//			normal_vector *= -1;
//	}
//
//	return normalized_normal_vectors;
//}
//

//std::vector<uint> Element::vertex_node_indexes(void) const {
	//const auto num_vertex = this->geometry_.reference_geometry_.num_vertex();

	//return { this->node_indexes_.begin(), this->node_indexes_.begin() + num_vertex };
//}
//

//
//


//
//

//

//Euclidean_Vector Element::node_at_index(const uint index) const {
//	const auto index_iter = std::find(this->node_indexes_.begin(), this->node_indexes_.end(), index);
//	dynamic_require(index_iter != this->node_indexes_.end(), "index should be included in node indexes");
//
//	const auto pos = index_iter - this->node_indexes_.begin();
//
//	const auto& this_nodes = this->geometry_.get_nodes();
//	return this_nodes[pos];
//}
//

//std::vector<Euclidean_Vector> Element::nodes_at_indexes(const std::vector<uint>& indexes) const {
//	const auto num_index = indexes.size();	
//	std::vector<Euclidean_Vector> nodes(num_index);
//
//	for (ushort i = 0; i < num_index; ++i) 
//		nodes[i] = this->node_at_index(indexes[i]);
//
//	return nodes;
//}
//



















namespace ms 
{
	template <typename T, typename Container>
	std::vector<T> extract_by_index(const std::vector<T>& set, const Container& indexes) 
	{
		const auto num_extracted_value = indexes.size();
		
		std::vector<T> extracted_values;
		extracted_values.reserve(num_extracted_value);

		for (const auto& index : indexes) 
			extracted_values.push_back(set[index]);

		return extracted_values;
	}

	template <typename T>
	bool is_circular_permutation(const std::vector<T>& set1, const std::vector<T>& set2) 
	{
		const auto set1_size = set1.size();
		const auto set2_size = set2.size();
				
		if (set1_size < set2_size)
			return false;

		std::vector<T> temp;
		temp.reserve(set1_size * 2);
		temp.insert(temp.end(), set1.begin(), set1.end());
		temp.insert(temp.end(), set1.begin(), set1.end());

		auto start_iter = temp.begin();
		for (size_t i = 0; i < set1_size; ++i, ++start_iter) 
		{
			if (std::equal(start_iter, start_iter + set2_size, set2.begin()))
				return true;
		}

		return false;
	}
//
//	template <typename T>
//	bool contains(const std::vector<T>& vec, const T& val) {
//		return std::find(vec.begin(), vec.end(), val) != vec.end();
//	}
//
//	template <typename T>
//	bool has_intersection(const std::vector<T>& vec1, const std::vector<T>& vec2) {
//		auto tmp1 = vec1;
//		auto tmp2 = vec2;
//
//		std::sort(tmp1.begin(), tmp1.end());
//		std::sort(tmp2.begin(), tmp2.end());
//
//		std::vector<T> intersection;
//		std::set_intersection(tmp1.begin(), tmp1.end(), tmp2.begin(), tmp2.end(), std::back_inserter(intersection));
//
//		return !intersection.empty();
//	}
}