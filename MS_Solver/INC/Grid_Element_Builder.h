#pragma once
#include "Grid_File_Type.h"
#include "Profiler.h"
#include "Log.h"

#include <map>
#include <unordered_set>
#include <sstream>


struct Grid_Elements
{
	std::vector<Element> cell_elements;
	std::vector<Element> boundary_elements;
	std::vector<std::pair<Element, Element>> periodic_boundary_element_pairs;
	std::vector<Element> inner_face_elements;
};

//
//template <ushort space_dimension>
//struct Grid_Elements
//{
	//std::vector<Element<space_dimension>> cell_elements;
	//std::vector<Element<space_dimension>> boundary_elements;
	//std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> periodic_boundary_element_pairs;
	//std::vector<Element<space_dimension>> inner_face_elements;
//};
//
//
//template <typename Grid_File_Type, ushort space_dimension>
//class Grid_Element_Builder;
//find_matched_periodic_node_indexes
//
//template <ushort space_dimension>
//class Grid_Element_Builder<Gmsh, space_dimension>
//{
//private:
//	Grid_Element_Builder(void) = delete;
//
//private:
//	using This_			= Grid_Element_Builder<Gmsh, space_dimension>;
//	using Euclidean_Vector = Euclidean_Vector<space_dimension>;
//
//public:
//	static Grid_Elements<space_dimension> build_from_grid_file(const std::string& grid_file_name);
//
//	//private: for test
//	static std::vector<Element<space_dimension>> make_inner_face_elements(const std::vector<Element<space_dimension>>& cell_elements, const std::vector<Element<space_dimension>>& boundary_elements, const std::vector<Element<space_dimension>>& periodic_boundary_elements);
//	static std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> match_periodic_boundaries(std::vector<Element<space_dimension>>& periodic_boundary_elements);
//};
//
//

//
//
////Template Definition
//template <ushort space_dimension>
//Grid_Elements<space_dimension> Grid_Element_Builder<Gmsh, space_dimension>::build_from_grid_file(const std::string& grid_file_name) {
	//SET_TIME_POINT;
	//Log::content_ << "================================================================================\n";
	//Log::content_ << "\t\t\t\t PreProcessing \n";
	//Log::content_ << "================================================================================\n";	
	//Log::print(); 		

	//SET_TIME_POINT;
	//
	//const auto grid_file_path = "RSC/Grid/" + std::to_string(space_dimension) + "D/" + grid_file_name + ".msh";


	//std::ifstream grid_file_stream(grid_file_path);
	//dynamic_require(grid_file_stream.is_open(), "fail to open " + grid_file_path);
	//
	//const auto node_text			= This_::read_about(grid_file_stream, "Nodes");
	//const auto node_datas			= This_::make_node_datas(node_text);
	//
	//const auto element_text			= This_::read_about(grid_file_stream, "Elements");	
	//const auto physical_name_text	= This_::read_about(grid_file_stream, "PhysicalNames");

	//Log::content_ << std::left << std::setw(50) << "@ Read Grid" << " ----------- " << GET_TIME_DURATION << "s\n\n";
	//Log::print();

	//return make_elements(element_text, physical_name_text, node_datas);
//}

//template <ushort space_dimension>
//Grid_Elements<space_dimension> Grid_Element_Builder<Gmsh, space_dimension>::make_elements(const Text& element_text, const Text& physical_name_text, const std::vector<Euclidean_Vector>& node_datas) {
//	SET_TIME_POINT;
//
//	std::map<ushort, ElementType> physical_group_index_to_element_type;
//
//	for (const auto& physical_name_sentence : physical_name_text) {
//		const char delimiter = ' ';
//		const auto parsed_sentence_set = ms::parse(physical_name_sentence, delimiter);
//
//		//const size_t dimension		= parsed_sentence_set[0].toValue<size_t>();
//		const auto physical_group_index	= ms::string_to_value<ushort>(parsed_sentence_set[1]);
//		const auto name					= ms::remove(parsed_sentence_set[2], "\"");
//		const auto element_type			= ms::string_to_element_type(name);
//
//		physical_group_index_to_element_type.emplace(physical_group_index, element_type);
//	}
//
//	std::vector<Element<space_dimension>> cell_elements;
//	std::vector<Element<space_dimension>> boundary_elements;
//	std::vector<Element<space_dimension>> periodic_boundary_elements;
//
//	for (const auto& element_sentence : element_text) {
//		const auto delimiter = ' ';
//		const auto parsed_sentences = ms::parse(element_sentence, delimiter);
//
//		auto value_set = ms::string_to_value_set<uint>(parsed_sentences);
//
//		//const auto index					= value_set[0];
//		const auto figure_type_index		= value_set[1];
//		//const auto tag_index				= value_set[2];
//		const auto physical_gorup_index		= value_set[3];
//		//const auto element_group_index	= value_set[4];
//
//		//reference geometry
//		const auto figure		= Gmsh::figure_type_index_to_element_figure(figure_type_index);
//		const auto figure_order = Gmsh::figure_type_index_to_figure_order(figure_type_index);
//		auto reference_geometry = ReferenceGeometry<space_dimension>(figure, figure_order);
//
//		//geometry
//		constexpr ushort num_index = 5;
//		value_set.erase(value_set.begin(), value_set.begin() + num_index);
//
//		const auto num_nodes = value_set.size();
//		std::vector<uint> node_indexes(num_nodes);
//		for (ushort i = 0; i < num_nodes; ++i)
//			node_indexes[i] = value_set[i] - 1;		//Gmsh node index start with 1
//
//		auto nodes = ms::extract_by_index(node_datas, node_indexes);
//		Geometry geometry(reference_geometry, std::move(nodes));
//
//		//element
//		const auto type	= physical_group_index_to_element_type.at(physical_gorup_index);
//		Element element(type, std::move(geometry), std::move(node_indexes));
//		
//		switch (type) {
//		case ElementType::cell:
//			cell_elements.emplace_back(std::move(element));
//			break;
//		case ElementType::periodic:
//			periodic_boundary_elements.emplace_back(std::move(element));
//			break;
//		default:
//			boundary_elements.emplace_back(std::move(element));
//			break;
//		}
//	}
//
//	const auto inner_face_elements = This_::make_inner_face_elements(cell_elements, boundary_elements, periodic_boundary_elements);
//	const auto periodic_boundary_element_pairs = This_::match_periodic_boundaries(periodic_boundary_elements);
//
//	Log::content_ << std::left << std::setw(50) << "@ Make Elements" << " ----------- " << GET_TIME_DURATION << "s\n";
//	Log::content_ << "  " << std::setw(8) << cell_elements.size() << " cell \n";
//	Log::content_ << "  " << std::setw(8) << boundary_elements.size() << " boundary\n";
//	Log::content_ << "  " << std::setw(8) << periodic_boundary_element_pairs.size() << " periodic boundary pair\n";
//	Log::content_ << "  " << std::setw(8) << inner_face_elements.size() << " inner face\n\n";
//	Log::print();
//	return { cell_elements, boundary_elements, periodic_boundary_element_pairs, inner_face_elements };
//}
//
//template<ushort space_dimension>
//std::vector<Element<space_dimension>> Grid_Element_Builder<Gmsh, space_dimension>::make_inner_face_elements(const std::vector<Element<space_dimension>>& cell_elements, const std::vector<Element<space_dimension>>& boundary_elements, const std::vector<Element<space_dimension>>& periodic_boundary_elements) {
//	//check constructed face
//	std::set<std::vector<uint>> constructed_face_vnode_index_set;
// 
//	for (const auto& boundray_element : boundary_elements) {
//		auto vnode_indexes = boundray_element.vertex_node_indexes();
//		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
//		constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
//	}
//	for (const auto& periodic_boundray_element : periodic_boundary_elements) {
//		auto vnode_indexes = periodic_boundray_element.vertex_node_indexes();
//		std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
//		constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
//	}
//
//	//construct face & select inner face
//	std::vector<Element<space_dimension>> inner_face_elements;
//
//	for (const auto& cell_element : cell_elements) {
//		auto face_elements = cell_element.make_face_elements();
//		for (auto& face_element : face_elements) {
//			auto vnode_indexes = face_element.vertex_node_indexes();
//			std::sort(vnode_indexes.begin(), vnode_indexes.end());	//to ignore index order
//
//			if (constructed_face_vnode_index_set.contains(vnode_indexes))
//				continue;
//			else {
//				constructed_face_vnode_index_set.insert(std::move(vnode_indexes));
//				inner_face_elements.push_back(std::move(face_element));
//			}
//		}
//	}
//
//	return inner_face_elements;
//}
//
//template<ushort space_dimension>
//std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> Grid_Element_Builder<Gmsh, space_dimension>::match_periodic_boundaries(std::vector<Element<space_dimension>>& periodic_boundary_elements) {	
//	const auto num_periodic_element = periodic_boundary_elements.size();
//	const auto num_pair = static_cast<uint>(0.5 * num_periodic_element);
//
//	std::unordered_set<uint> matched_index_set;
//	matched_index_set.reserve(num_periodic_element);
//
//	std::vector<std::pair<Element<space_dimension>, Element<space_dimension>>> matched_periodic_element_pairs;
//	matched_periodic_element_pairs.reserve(num_pair);
//
//	for (uint i = 0; i < num_periodic_element; ++i) {
//		if (matched_index_set.contains(i))
//			continue;
//
//		auto& i_element = periodic_boundary_elements[i];
//
//		for (uint j = i + 1; j < num_periodic_element; ++j) {
//			if (matched_index_set.contains(j))
//				continue;
//			
//			const auto& j_element = periodic_boundary_elements[j];
//
//			auto matched_periodic_node_indexes = i_element.find_matched_periodic_node_indexes(j_element);
//
//			if (!matched_periodic_node_indexes.empty()) {
//				//modify j_element based on i_element				
//				auto matched_periodic_nodes = j_element.nodes_at_indexes(matched_periodic_node_indexes);
//
//				const auto j_reference_geometry = j_element.geometry_.reference_geometry_;
//				Geometry matched_j_geometry(j_reference_geometry, std::move(matched_periodic_nodes));
//				Element matched_j_element(ElementType::periodic, std::move(matched_j_geometry), std::move(matched_periodic_node_indexes));
//
//				matched_periodic_element_pairs.push_back(std::make_pair(std::move(i_element), std::move(matched_j_element)));
//				matched_index_set.insert(i);
//				matched_index_set.insert(j);
//				break;
//			}
//		}
//	}
//
//	dynamic_require(matched_periodic_element_pairs.size() == num_pair, "every pbdry should be matched");
//	return matched_periodic_element_pairs;
//}



////inline function definition
//namespace ms {
//	inline ElementType string_to_element_type(const std::string& str) {
//		if (ms::contains_icase(str, "Unspecified"))
//			return ElementType::cell;
//		else if (ms::contains_icase(str, "slip", "wall"))
//			return ElementType::slip_wall;
//		else if (ms::contains_icase(str, "SuperSonic", "inlet", "1"))
//			return ElementType::supersonic_inlet1;
//		else if (ms::contains_icase(str, "SuperSonic", "inlet", "2"))
//			return ElementType::supersonic_inlet2;
//		else if (ms::contains_icase(str, "SuperSonic", "outlet"))
//			return ElementType::supersonic_outlet;
//		else if (ms::contains_icase(str, "periodic"))
//			return ElementType::periodic;
//		else {
//			throw std::runtime_error("wrong element_type");
//			return ElementType::not_in_list;
//		}
//	}
//}
