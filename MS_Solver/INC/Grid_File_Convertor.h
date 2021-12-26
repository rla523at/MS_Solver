#pragma once
#include "Configuration.h"
#include "Element.h"

#include "Log.h"
#include "Profiler.h"

class Grid_File_Convertor
{
public:
	virtual std::vector<Element> convert_to_elements(const std::string_view grid_file_name) const abstract;
	ushort get_space_dimension(void) const
	{
		return this->space_dimension_;
	};

protected:
	ushort space_dimension_ = 0;
};

class Gmsh_Convertor : public Grid_File_Convertor
{
	enum class GmshFigureType
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

public:
	Gmsh_Convertor(const ushort space_dimension);

public:
	std::vector<Element> convert_to_elements(const std::string_view grid_file_path) const;

private:
	Text read_about(std::ifstream& grid_file_stream, const std::string& target) const;
	std::vector<Euclidean_Vector> make_nodes(const Text& node_text) const;		
	std::map<ushort, ElementType> make_physical_group_index_to_element_type(const Text& physical_name_text) const;
	std::vector<Element> make_elements(const Text& element_text, const std::map<ushort, ElementType>& physical_group_index_to_element_type, const std::vector<Euclidean_Vector>& node_datas) const;
	ushort figure_type_index_to_figure_order(const ushort figure_type_index) const;
	Figure figure_type_index_to_element_figure(const ushort figure_type_index) const;
};

class Grid_File_Convertor_Factory//static class
{
public:
	static std::unique_ptr<Grid_File_Convertor> make_unique(const Configuration& configuration);

private:
	Grid_File_Convertor_Factory(void) = delete;
};


namespace ms 
{
	ElementType sentece_to_element_type(const Sentence& sentence);
}