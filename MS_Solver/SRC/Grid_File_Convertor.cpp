#include "../INC/Grid_File_Convertor.h"

Gmsh_Convertor::Gmsh_Convertor(const ushort space_dimension)
{
	this->space_dimension_ = space_dimension;
}


std::vector<Element> Gmsh_Convertor::convert_to_elements(const std::string_view grid_file_path) const
{
	Profiler::set_time_point();

	std::ifstream grid_file_stream(grid_file_path);
	REQUIRE(grid_file_stream.is_open(), "fail to open grid file");

	const auto node_text = this->read_about(grid_file_stream, "Nodes");
	const auto element_text = this->read_about(grid_file_stream, "Elements");
	const auto physical_name_text = this->read_about(grid_file_stream, "PhysicalNames");

	LOG << std::left << std::setw(50) << "@ Read Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;

	Profiler::set_time_point();

	const auto nodes = this->make_nodes(node_text);
	const auto physical_group_index_to_element_type = this->make_physical_group_index_to_element_type(physical_name_text);
	return this->make_elements(element_text, physical_group_index_to_element_type, nodes);


	//복사가 발생해서 error가 나는데 왜 어디서 복사가 발생하는거지 ?
	//const auto elements = this->make_elements(element_text, physical_group_index_to_element_type, nodes);
	//return elements;
	//LOG << std::left << std::setw(50) << "@ Convert Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;

}

Text Gmsh_Convertor::read_about(std::ifstream& grid_file_stream, const std::string& target) const
{
	bool is_find = false;
	const auto target_str = "$" + target;

	std::string tmp_str;
	while (std::getline(grid_file_stream, tmp_str)) 
	{		
		if (ms::contains_icase(tmp_str,target_str)) 
		{
			is_find = true;
			std::getline(grid_file_stream, tmp_str);
			break;
		}
	}

	REQUIRE(is_find, "\"" + target_str + "\" should be contained in grid file");

	const auto num_data = ms::string_to_value<size_t>(tmp_str);

	Text txt;
	txt.read(grid_file_stream, num_data);
	return txt;
}

std::vector<Euclidean_Vector> Gmsh_Convertor::make_nodes(const Text& node_text) const
{
	std::vector<Euclidean_Vector> node_datas;
	node_datas.reserve(node_text.size());

	for (const auto& node_sentence : node_text)
	{
		const char delimiter = ' ';
		auto parsed_sentences = node_sentence.parse(delimiter);

		std::vector<double> node_coords(this->space_dimension_);
		for (ushort i = 0; i < this->space_dimension_; ++i)
		{
			node_coords[i] = parsed_sentences[i + 1].to_value<double>();// parsed_sentences[0] : node index
		}

		node_datas.push_back(std::move(node_coords));
	}

	return node_datas;
}

std::map<ushort, ElementType> Gmsh_Convertor::make_physical_group_index_to_element_type(const Text& physical_name_text) const
{
	std::map<ushort, ElementType> physical_group_index_to_element_type;

	for (const auto& physical_name_sentence : physical_name_text)
	{
		const char delimiter = ' ';
		const auto parsed_sentences = physical_name_sentence.parse(delimiter);

		//const size_t dimension		= parsed_sentence_set[0].to_value<size_t>();
		const auto physical_group_index = parsed_sentences[1].to_value<ushort>();
		const auto name					= parsed_sentences[2].get_remove("\"");
		const auto element_type			= ms::sentece_to_element_type(name);

		physical_group_index_to_element_type.emplace(physical_group_index, element_type);
	}

	return physical_group_index_to_element_type;
}


std::vector<Element> Gmsh_Convertor::make_elements(const Text& element_text, const std::map<ushort, ElementType>& physical_group_index_to_element_type, const std::vector<Euclidean_Vector>& node_datas) const
{
	std::vector<Element> elements;
	elements.reserve(element_text.size());

	for (const auto& element_sentence : element_text)
	{
		const auto delimiter = ' ';
		const auto parsed_sentences = element_sentence.parse(delimiter);

		auto value_set = ms::sentences_to_values<uint>(parsed_sentences);

		//const auto index					= value_set[0];
		const auto figure_type_index		= value_set[1];
		//const auto tag_index				= value_set[2];
		const auto physical_gorup_index		= value_set[3];
		//const auto element_group_index	= value_set[4];

		//reference geometry
		const auto figure = this->figure_type_index_to_element_figure(figure_type_index);
		const auto figure_order = this->figure_type_index_to_figure_order(figure_type_index);
		auto reference_geometry = Reference_Geometry_Factory::make_unique(figure, figure_order);

		//geometry
		constexpr auto num_indexes = 5;
		value_set.erase(value_set.begin(), value_set.begin() + num_indexes);

		const auto num_nodes = value_set.size();
		std::vector<uint> node_indexes(num_nodes);
		for (ushort i = 0; i < num_nodes; ++i)
		{
			node_indexes[i] = value_set[i] - 1;		//Gmsh node index start with 1
		}

		auto nodes = ms::extract_by_index(node_datas, node_indexes);
		Geometry geometry(std::move(reference_geometry), std::move(nodes));

		//element
		const auto type = physical_group_index_to_element_type.at(physical_gorup_index);
		elements.push_back({ type, std::move(node_indexes), std::move(geometry) });
	}
	
	LOG << std::left << std::setw(50) << "@ Convert Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;

	return elements;
}

ushort Gmsh_Convertor::figure_type_index_to_figure_order(const ushort element_type_indx) const
{
	switch (static_cast<GmshFigureType>(element_type_indx)) 
	{
	case GmshFigureType::POINT:			return 0;
	case GmshFigureType::LINE_P1:
	case GmshFigureType::TRIS_P1:
	case GmshFigureType::QUAD_P1:
	case GmshFigureType::TETS_P1:
	case GmshFigureType::HEXA_P1:
	case GmshFigureType::PRIS_P1:
	case GmshFigureType::PYRA_P1:		return 1;
	case GmshFigureType::LINE_P2:
	case GmshFigureType::TRIS_P2:
	case GmshFigureType::QUAD_P2:
	case GmshFigureType::TETS_P2:
	case GmshFigureType::HEXA_P2:
	case GmshFigureType::PRIS_P2:
	case GmshFigureType::PYRA_P2:		return 2;
	case GmshFigureType::LINE_P3:
	case GmshFigureType::TRIS_P3:
	case GmshFigureType::QUAD_P3:
	case GmshFigureType::TETS_P3:
	case GmshFigureType::HEXA_P3:
	case GmshFigureType::PRIS_P3:
	case GmshFigureType::PYRA_P3:		return 3;
	case GmshFigureType::LINE_P4:
	case GmshFigureType::TRIS_P4:
	case GmshFigureType::QUAD_P4:
	case GmshFigureType::TETS_P4:
	case GmshFigureType::HEXA_P4:
	case GmshFigureType::PRIS_P4:
	case GmshFigureType::PYRA_P4:		return 4;
	case GmshFigureType::LINE_P5:
	case GmshFigureType::TRIS_P5:
	case GmshFigureType::QUAD_P5:
	case GmshFigureType::TETS_P5:
	case GmshFigureType::HEXA_P5:
	case GmshFigureType::PRIS_P5:		return 5;
	case GmshFigureType::LINE_P6:
	case GmshFigureType::QUAD_P6:		return 6;
	default:
		EXCEPTION("invalid element type index");
		return NULL;
	}
}


Figure Gmsh_Convertor::figure_type_index_to_element_figure(const ushort element_type_index) const
{
	switch (static_cast<GmshFigureType>(element_type_index)) 
	{
	case GmshFigureType::POINT:			return Figure::point;
	case GmshFigureType::LINE_P1:
	case GmshFigureType::LINE_P2:
	case GmshFigureType::LINE_P3:
	case GmshFigureType::LINE_P4:
	case GmshFigureType::LINE_P5:
	case GmshFigureType::LINE_P6:		return Figure::line;
	case GmshFigureType::TRIS_P1:
	case GmshFigureType::TRIS_P2:
	case GmshFigureType::TRIS_P3:
	case GmshFigureType::TRIS_P4:
	case GmshFigureType::TRIS_P5:		return Figure::triangle;
	case GmshFigureType::QUAD_P1:
	case GmshFigureType::QUAD_P2:
	case GmshFigureType::QUAD_P3:
	case GmshFigureType::QUAD_P4:
	case GmshFigureType::QUAD_P5:
	case GmshFigureType::QUAD_P6:		return Figure::quadrilateral;
	case GmshFigureType::TETS_P1:
	case GmshFigureType::TETS_P2:
	case GmshFigureType::TETS_P3:
	case GmshFigureType::TETS_P4:
	case GmshFigureType::TETS_P5:		return Figure::tetrahedral;
	case GmshFigureType::HEXA_P1:
	case GmshFigureType::HEXA_P2:
	case GmshFigureType::HEXA_P3:
	case GmshFigureType::HEXA_P4:
	case GmshFigureType::HEXA_P5:		return Figure::hexahedral;
	case GmshFigureType::PRIS_P1:
	case GmshFigureType::PRIS_P2:
	case GmshFigureType::PRIS_P3:
	case GmshFigureType::PRIS_P4:
	case GmshFigureType::PRIS_P5:		return Figure::prism;
	case GmshFigureType::PYRA_P1:
	case GmshFigureType::PYRA_P2:
	case GmshFigureType::PYRA_P3:
	case GmshFigureType::PYRA_P4:		return Figure::pyramid;
	default:
		EXCEPTION("invalid element type index");
		return Figure::point;
	}
}

std::unique_ptr<Grid_File_Convertor> Grid_File_Convertor_Factory::make_unique(const Configuration& configuration)
{
	const auto space_dimension = configuration.get<ushort>("space_dimension");
	const auto grid_file_type = configuration.get("grid_file_type");

	if (ms::contains_icase(grid_file_type, "gmsh"))
	{
		return std::make_unique<Gmsh_Convertor>(space_dimension);
	}
	else
	{
		EXCEPTION("grid file type in configuration file does not supported");
		return nullptr;
	}
}

namespace ms {
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