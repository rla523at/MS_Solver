#pragma once
#include "Figure.h"
#include "Grid_Data.h"
#include "Text.h"

#include "Log.h"
#include "Profiler.h"

#include <map>
#include <random>
#include <unordered_set>

using ushort = unsigned short;

class Gmsh_Grid_File
{
private:
	enum class Gmsh_Figure_Type
	{
		POINT = 0,
		LINE_P1 = 1, LINE_P2 = 8, LINE_P3 = 26, LINE_P4 = 27, LINE_P5 = 28, LINE_P6 = 62,
		TRIS_P1 = 2, TRIS_P2 = 9, TRIS_P3 = 21, TRIS_P4 = 23, TRIS_P5 = 25,
		QUAD_P1 = 3, QUAD_P2 = 10, QUAD_P3 = 36, QUAD_P4 = 37, QUAD_P5 = 38, QUAD_P6 = 47,
		TETS_P1 = 4, TETS_P2 = 11, TETS_P3 = 29, TETS_P4 = 30, TETS_P5 = 31,
		HEXA_P1 = 5, HEXA_P2 = 12, HEXA_P3 = 92, HEXA_P4 = 93, HEXA_P5 = 94,
		PRIS_P1 = 6, PRIS_P2 = 13, PRIS_P3 = 90, PRIS_P4 = 91, PRIS_P5 = 106,
		PYRA_P1 = 7, PYRA_P2 = 14, PYRA_P3 = 118, PYRA_P4 = 119
	};
	enum class Gmsh_Data_List
	{
		node_data, element_data, physical_name_data
	};

public:
	Gmsh_Grid_File(const short space_dimension, const std::string_view grid_file_path);

public://Command
	void perturb_node(const double maximum_perturbation_size);

public://Query 
	std::vector<Euclidean_Vector> extract_node_datas(void) const;
	std::vector<Element_Data> extract_element_datas(void) const;	
	void write(const std::string_view output_file_path) const;

private:
	ElementType convert_to_element_type(const Sentence& name) const;
	Figure convert_to_figure(const ushort figure_type_index) const;
	short convert_to_figure_order(const ushort figure_type_index) const;
	Text extract_text_about(const Gmsh_Data_List target_data) const;
	std::pair<size_t, size_t> find_start_end_line_index_pair(const Gmsh_Data_List target_data) const;
	std::unordered_set<size_t> find_bdry_node_index_set(void) const;
	std::string gmsh_data_to_keyword(const Gmsh_Data_List target_data) const;
	std::map<ushort, ElementType> make_physical_group_index_to_element_type(void) const;

private:
	ushort space_dimension_;
	Text grid_text_;
};