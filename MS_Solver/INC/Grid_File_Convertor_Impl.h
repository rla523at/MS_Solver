//#pragma once
//#include "Grid_File_Convertor.h"
//#include "Grid_File_Reader.h"
//
//#include "Configuration.h"
//#include "Reference_Geometry_Impl.h"
//#include "Log.h"
//#include "Profiler.h"
//
//
//class Gmsh_Convertor : public Grid_File_Convertor
//{
//public:
//	Gmsh_Convertor(const ushort space_dimension);
//
//public:
//	std::vector<Element> convert_to_elements(const std::string_view grid_file_path) const;
//
//private:
//	std::vector<Euclidean_Vector> make_nodes(const Text& node_text) const;
//	std::map<ushort, ElementType> make_physical_group_index_to_element_type(const Text& physical_name_text) const;
//	std::vector<Element> make_elements(const Text& element_text, const std::map<ushort, ElementType>& physical_group_index_to_element_type, const std::vector<Euclidean_Vector>& node_datas) const;
//
//private:
//	Gmsh_File_Reader gmsh_file_reader_;
//};
//
//class Grid_File_Convertor_Factory//static class
//{
//public:
//	static std::unique_ptr<Grid_File_Convertor> make_unique(const Configuration& configuration);
//
//private:
//	Grid_File_Convertor_Factory(void) = delete;
//};