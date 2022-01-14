#pragma once
#include "Figure.h"
#include "Grid_File_Type.h"
#include "Grid_Data.h"
#include "Text.h"

#include "Log.h"
#include "Profiler.h"

namespace ms
{
	ElementType sentece_to_element_type(const Sentence& sentence)
	{
		if (sentence.contain_icase("Unspecified"))
		{
			return ElementType::cell;
		}
		else if (sentence.contain_icase("slip", "wall"))
		{
			return ElementType::slip_wall;
		}
		else if (sentence.contain_icase("SuperSonic", "inlet", "1"))
		{
			return ElementType::supersonic_inlet1;
		}
		else if (sentence.contain_icase("SuperSonic", "inlet", "2"))
		{
			return ElementType::supersonic_inlet2;
		}
		else if (sentence.contain_icase("SuperSonic", "outlet"))
		{
			return ElementType::supersonic_outlet;
		}
		else if (sentence.contain_icase("periodic"))
		{
			return ElementType::periodic;
		}
		else if (sentence.contain_icase("initial", "constant"))
		{
			return ElementType::initial_constant_BC;
		}
		else
		{
			EXCEPTION("wrong element_type");
			return ElementType::not_in_list;
		}
	};
}

class Gmsh_File_Reader
{
public:
	Gmsh_File_Reader(const short space_dimension)
	{
		REQUIRE(0 < space_dimension, "space dimension should be positive ");
		REQUIRE(space_dimension <= 3, "space dimension can not exceed 3");

		this->space_dimension_ = space_dimension;
	}

public:
	std::vector<Euclidean_Vector> read_node_datas(const Text& grid_text) const
	{
		constexpr auto keyword = "$Nodes";
		const auto node_text = this->extract_about(grid_text, keyword);

		const auto num_node_datas = node_text.size(); 
		std::vector<Euclidean_Vector> node_datas(num_node_datas + 1); //To match node index and node data index (* Gmsh node index start with 1)

		for (size_t i = 0; i < num_node_datas; ++i)
		{
			auto& data = node_datas[i + 1]; //Gmsh node index start with 1

			const auto& sentece = node_text[i];

			const char delimiter = ' ';
			auto parsed_sentences = sentece.parse(delimiter);

			data = Euclidean_Vector(this->space_dimension_);
			for (ushort i = 0; i < this->space_dimension_; ++i)
			{
				data[i] = parsed_sentences[i + 1].to_value<double>();// parsed_sentences[0] : node index
			}
		}

		return node_datas;
	}

	std::vector<Element_Data> read_element_datas(const Text& grid_text) const
	{
		const auto physical_group_index_to_element_type = this->make_physical_group_index_to_element_type(grid_text);

		constexpr auto keyword = "$Elements";
		const auto element_text = this->extract_about(grid_text, keyword);

		const auto num_elements = element_text.size();
		
		std::vector<Element_Data> element_datas(num_elements);

		for (size_t i = 0; i < num_elements; ++i)
		{
			auto& data = element_datas[i];
			
			const auto& element_sentence = element_text[i];

			const auto delimiter = ' ';
			const auto parsed_sentences = element_sentence.parse(delimiter);

			auto value_set = ms::sentences_to_values<uint>(parsed_sentences);

			//const auto index					= value_set[0];
			const auto figure_type_index = value_set[1];
			//const auto tag_index				= value_set[2];
			const auto physical_gorup_index = value_set[3];
			//const auto element_group_index	= value_set[4];

			data.figure = Gmsh::figure_type_index_to_figure(figure_type_index);
			data.figure_order = Gmsh::figure_type_index_to_figure_order(figure_type_index);
			data.element_type = physical_group_index_to_element_type.at(physical_gorup_index);

			constexpr auto num_indexes = 5;
			value_set.erase(value_set.begin(), value_set.begin() + num_indexes);

			const auto num_nodes = value_set.size();
			data.node_indexes.resize(num_nodes);
			
			for (ushort i = 0; i < num_nodes; ++i)
			{
				data.node_indexes[i] = value_set[i];		
			}
		}
		return element_datas;
	}

	std::pair<size_t, size_t> node_data_start_end_line_index_pair(const Text& grid_text) const
	{
		return this->data_start_end_line_index_pair(grid_text, "$Nodes");
	}

private:
	std::pair<size_t, size_t> data_start_end_line_index_pair(const Text& grid_text, const std::string_view keyword) const
	{
		const auto line_index = grid_text.find_line_index_has_keyword(keyword);

		const auto& sentece = grid_text[line_index + 1];
		const auto num_data = sentece.to_value<size_t>();

		const auto start_line_index = line_index + 2;
		const auto end_line_index = start_line_index + num_data;

		return { start_line_index, end_line_index };
	}

	Text extract_about(const Text& grid_text, const std::string_view keyword) const
	{
		const auto [start_line_index, end_line_index] = this->data_start_end_line_index_pair(grid_text, keyword);
		return grid_text.extract(start_line_index, end_line_index);
	}

	std::map<ushort, ElementType> make_physical_group_index_to_element_type(const Text& grid_text) const
	{
		constexpr auto keyword = "$PhysicalNames";
		const auto physical_name_text = this->extract_about(grid_text, keyword);

		std::map<ushort, ElementType> physical_group_index_to_element_type;

		for (const auto& physical_name_sentence : physical_name_text)
		{
			const char delimiter = ' ';
			const auto parsed_sentences = physical_name_sentence.parse(delimiter);

			//const size_t dimension		= parsed_sentence_set[0].to_value<size_t>();
			const auto physical_group_index = parsed_sentences[1].to_value<ushort>();
			const auto name = parsed_sentences[2].get_remove("\"");
			const auto element_type = ms::sentece_to_element_type(name);

			physical_group_index_to_element_type.emplace(physical_group_index, element_type);
		}

		return physical_group_index_to_element_type;
	}

private:
	ushort space_dimension_;		
};