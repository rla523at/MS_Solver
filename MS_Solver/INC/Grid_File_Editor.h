#pragma once
#include "Grid_File_Reader.h"
#include "Element_Type.h"

#include <unordered_set>
#include <random>

class Gmsh_Text_Editor
{
public:
	Gmsh_Text_Editor(const short space_dimension)
		: file_reader_(space_dimension)
		, space_dimension_(space_dimension) {};

public:
	void perturb(Text& grid_text, const double characteristic_length, const ushort max_pertubation_percentage)
	{		
		REQUIRE(max_pertubation_percentage <= 100, "percentage can not exceed 100%");

		const auto scail_factor = static_cast<double>(max_pertubation_percentage) / 100.0;
		const auto perturbation_size = characteristic_length * scail_factor;

		std::random_device rd;
		std::mt19937 gen(rd());
		std::uniform_real_distribution<double> dis(-perturbation_size, perturbation_size);

		const auto bdry_node_index_set = this->find_bdry_node_index_set(grid_text);

		const auto [start_line_index, end_line_index] = this->file_reader_.node_data_start_end_line_index_pair(grid_text);
		
		for (size_t i = start_line_index; i < end_line_index; ++i)
		{
			auto& node_sentence = grid_text[i];

			constexpr auto denominator = ' ';
			const auto parsed_sentence = node_sentence.parse(denominator);

			const auto node_index = parsed_sentence[0].to_value<size_t>();

			if (bdry_node_index_set.contains(node_index))
			{
				continue;
			}
			else
			{
				constexpr auto num_coordinate = 3;

				node_sentence.clear();
				node_sentence << node_index;

				for (ushort j = 0; j < this->space_dimension_; ++j)
				{
					const auto coordinate = parsed_sentence[1 + j].to_value<double>(); // [0] is index
					const auto perturbation = dis(gen);
					const auto perturbed_coordinate = coordinate + perturbation;

					node_sentence.insert_with_space(perturbed_coordinate);
				}

				for (ushort j = 0; j < num_coordinate - this->space_dimension_; ++j)
				{
					constexpr auto zero = 0;
					node_sentence.insert_with_space(zero);
				}

			}
		}
	}

private:
	std::unordered_set<size_t> find_bdry_node_index_set(const Text& grid_text) const
	{
		const auto element_datas = this->file_reader_.read_element_datas(grid_text);
		
		const auto num_element = element_datas.size();

		std::unordered_set<size_t> bdry_node_index_set;

		for (size_t i = 0; i < num_element; ++i)
		{
			const auto& data = element_datas[i];

			if (data.element_type != ElementType::cell)
			{				
				bdry_node_index_set.insert(data.node_indexes.begin(), data.node_indexes.end());
			}
		}

		return bdry_node_index_set;
	}

private:
	ushort space_dimension_;
	Gmsh_File_Reader file_reader_;
};