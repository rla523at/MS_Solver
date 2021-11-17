#include "../INC/Element.h"

void Element::rearrange_node_indexes(std::vector<uint>&& rearranged_node_indexes)
{
	auto nodes = this->nodes_at_indexes(rearranged_node_indexes);
	this->node_indexes_ = std::move(rearranged_node_indexes);
	this->change_nodes(std::move(nodes));
}

bool Element::operator==(const Element& other) const
{
	return this->element_type_ == other.element_type_ &&
		this->node_indexes_ == other.node_indexes_ &&
		*this->reference_geometry_ == *other.reference_geometry_;
}


std::vector<uint> Element::find_periodic_matched_node_indexes(const Element& other) const 
{
	REQUIRE(this->element_type_ == ElementType::periodic && other.element_type_ == ElementType::periodic, "both element should be periodic");

	if (!this->can_be_periodic_pair(other))
		return {};

	const auto this_num_node = this->node_indexes_.size();
	const auto other_num_node = other.node_indexes_.size();

	std::unordered_set<uint> matched_vnode_index;
	matched_vnode_index.reserve(other_num_node);

	std::vector<uint> matched_periodic_node_indexes(this_num_node);

	const auto& this_nodes = this->nodes_;
	const auto& other_nodes = other.nodes_;

	for (int i = 0; i < other_num_node; ++i)
	{
		const auto& other_node = other_nodes[i];
		const auto other_node_index = other.node_indexes_[i];

		for (int j = 0; j < this_num_node; ++j) 
		{
			const auto& this_node = this_nodes[j];
			const auto this_node_index = this->node_indexes_[j];

			if (matched_vnode_index.contains(other_node_index))
			{
				continue;
			}

			if (this_node.is_axis_translation(other_node)) 
			{
				matched_periodic_node_indexes[i] = this_node_index;
				matched_vnode_index.insert(this_node_index);
				break;
			}
		}

		//when can not find pair
		if (i + 1 != matched_vnode_index.size())
			return {};
	}

	return matched_periodic_node_indexes;
}

bool Element::is_periodic_pair(const Element& other) const 
{
	REQUIRE(this->is_periodic_boundary() && other.is_periodic_boundary(), "both elemets should be periodic boundary");

	if (this->element_type_ != other.element_type_)
		return false;

	if (this->can_be_periodic_pair(other))
		return true;
	else
		return false;
}

std::vector<Element> Element::make_face_elements(void) const 
{
	REQUIRE(this->element_type_ == ElementType::cell, "make inner face elements should be called from cell element");

	auto face_geometries = this->face_geometries();
	auto set_of_face_node_indexes = this->set_of_face_node_indexes();

	const auto num_face = face_geometries.size();
	std::vector<Element> inner_face_elements;
	inner_face_elements.reserve(num_face);

	for (ushort i = 0; i < num_face; ++i)
		inner_face_elements.push_back({ ElementType::face, std::move(set_of_face_node_indexes[i]), std::move(face_geometries[i]) });

	return inner_face_elements;
}

std::vector<Euclidean_Vector> Element::nodes_at_indexes(const std::vector<uint>& indexes) const 
{
	const auto num_index = indexes.size();
	std::vector<Euclidean_Vector> nodes(num_index);

	for (ushort i = 0; i < num_index; ++i)
	{
		const auto index_iter = std::find(this->node_indexes_.begin(), this->node_indexes_.end(), indexes[i]);
		REQUIRE(index_iter != this->node_indexes_.end(), "index should be included");

		const auto pos = index_iter - this->node_indexes_.begin();
		nodes[i] = this->nodes_[pos];
	}
		
	return nodes;
}

Euclidean_Vector Element::outward_normalized_normal_vector(const Element& owner_cell_element, const Euclidean_Vector& node) const
{
	auto normal_vector = this->normalized_normal_vector(node);

	const auto face_type = this->check_face_type(owner_cell_element);	
	if (face_type == FaceType::inward_face)
	{
		normal_vector *= -1.0;
	}

	return normal_vector;
}

std::vector<Euclidean_Vector> Element::outward_normalized_normal_vectors(const Element& owner_cell_element, const std::vector<Euclidean_Vector>& points) const
{
	auto normalized_normal_vectors = this->normalized_normal_vectors(points);
	
	const auto face_type = this->check_face_type(owner_cell_element);
	if (face_type == FaceType::inward_face) 
	{
		for (auto& normal_vector : normalized_normal_vectors)
			normal_vector *= -1;
	}

	return normalized_normal_vectors;
}


ElementType Element::type(void) const 
{
	return this->element_type_;
}

std::vector<uint> Element::vertex_node_indexes(void) const
{
	const auto num_vertex = this->reference_geometry_->num_vertex();

	return { this->node_indexes_.begin(), this->node_indexes_.begin() + num_vertex };
}

bool Element::can_be_periodic_pair(const Element& other) const
{
	if (*this->reference_geometry_ != *other.reference_geometry_)
		return false;

	if (this->is_on_same_axis_plane(other))
		return false;

	return true;
}

FaceType Element::check_face_type(const Element& owner_cell_element) const 
{
	REQUIRE(this->element_type_ != ElementType::cell, "face or boundary element should be use this method");

	const auto this_vnode_indexes = this->vertex_node_indexes();
	const auto set_of_face_vnode_indexes = owner_cell_element.set_of_face_vertex_node_indexes();

	for (const auto& face_vnode_indexes : set_of_face_vnode_indexes) 
	{
		if (!std::is_permutation(face_vnode_indexes.begin(), face_vnode_indexes.end(), this_vnode_indexes.begin()))
			continue;

		if (this->is_line()) 
		{
			if (this_vnode_indexes == face_vnode_indexes)
				return FaceType::inward_face;
			else
				return FaceType::outward_face;
		}
		else
		{
			if (ms::is_circular_permutation(face_vnode_indexes, this_vnode_indexes))
				return FaceType::inward_face;
			else
				return FaceType::outward_face;
		}
	}

	EXCEPTION("Owner Cell : this face is not my face!");
	return FaceType::not_my_face;
}

bool Element::is_periodic_boundary(void) const 
{
	return this->element_type_ == ElementType::periodic;
}

std::vector<std::vector<uint>> Element::set_of_face_node_indexes(void) const 
{
	const auto set_of_face_node_index_sequences = this->reference_geometry_->set_of_face_node_index_sequences();
	const auto num_face = set_of_face_node_index_sequences.size();

	std::vector<std::vector<uint>> set_of_face_node_indexes(num_face);
	for (uint i = 0; i < num_face; ++i)
	{
		set_of_face_node_indexes[i] = ms::extract_by_index(this->node_indexes_, set_of_face_node_index_sequences[i]);
	}

	return set_of_face_node_indexes;
}

std::vector<std::vector<uint>> Element::set_of_face_vertex_node_indexes(void) const 
{
	const auto set_of_face_vnode_index_sequences = this->reference_geometry_->set_of_face_vertex_node_index_sequences();
	const auto num_face = set_of_face_vnode_index_sequences.size();

	std::vector<std::vector<uint>> set_of_face_vnode_indexes(num_face);
	for (uint i = 0; i < num_face; ++i)
	{
		set_of_face_vnode_indexes[i] = ms::extract_by_index(this->node_indexes_, set_of_face_vnode_index_sequences[i]);
	}

	return set_of_face_vnode_indexes;
}