#include "../INC/Grid_File.h"

Gmsh_Grid_File::Gmsh_Grid_File(const short space_dimension, const std::string_view grid_file_path)
{
	REQUIRE(0 < space_dimension, "space dimension should be positive ");
	REQUIRE(space_dimension <= 3, "space dimension can not exceed 3");

	Profiler::set_time_point();

	this->space_dimension_ = space_dimension;
	this->grid_text_.read(grid_file_path);

	LOG << std::left << std::setw(50) << "@ Read Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;
}

void Gmsh_Grid_File::perturb_node(const double maximum_perturbation_size)
{
	REQUIRE(max_pertubation_percentage <= 100, "percentage can not exceed 100%");

	std::random_device rd;
	std::mt19937 gen(rd());
	std::uniform_real_distribution<double> dis(-maximum_perturbation_size, maximum_perturbation_size);

	const auto bdry_node_index_set = this->find_bdry_node_index_set();

	const auto [start_line_index, end_line_index] = this->find_start_end_line_index_pair(Gmsh_Data_List::node_data);

	for (size_t i = start_line_index; i < end_line_index; ++i)
	{
		auto& node_sentence = this->grid_text_[i];

		constexpr auto denominator = ' ';
		const auto parsed_sentence = node_sentence.parse_by(denominator);

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
				const auto index = 1 + j;

				const auto coordinate = parsed_sentence[index].to_value<double>(); // [0] is node index
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

std::vector<Euclidean_Vector> Gmsh_Grid_File::extract_node_datas(void) const
{	
	Profiler::set_time_point();

	const auto node_data_text = this->extract_text_about(Gmsh_Data_List::node_data);

	const auto num_node_datas = node_data_text.size();
	std::vector<Euclidean_Vector> node_datas(num_node_datas + 1); //To match node index and node data index (* Gmsh node index start with 1)

	for (size_t i = 0; i < num_node_datas; ++i)
	{
		const auto& sentece = node_data_text[i];

		const char delimiter = ' ';
		const auto parsed_sentences = sentece.parse_by(delimiter);

		const auto node_index = parsed_sentences[0].to_value<size_t>();
		auto& data = node_datas[node_index];
		data = Euclidean_Vector(this->space_dimension_);

		for (ushort j = 0; j < this->space_dimension_; ++j)
		{
			data[j] = parsed_sentences[j + 1].to_value<double>();// parsed_sentences[0] : node index
		}
	}

	LOG << std::left << std::setw(50) << "@ Extract " << num_node_datas <<  " Node Data " << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;

	return node_datas;
}

std::vector<Element_Data> Gmsh_Grid_File::extract_element_datas(void) const
{
	Profiler::set_time_point();

	const auto physical_group_index_to_element_type = this->make_physical_group_index_to_element_type();

	const auto element_text = this->extract_text_about(Gmsh_Data_List::element_data);
	const auto num_elements = element_text.size();

	std::vector<Element_Data> element_datas(num_elements);

	for (size_t i = 0; i < num_elements; ++i)
	{
		const auto& sentence = element_text[i];
		const auto parsed_sentences = sentence.parse_by(' ');

		auto indexes = ms::sentences_to_values<uint>(parsed_sentences);

		//const auto index					= value_set[0];
		const auto figure_type_index = indexes[1];
		//const auto tag_index				= value_set[2];
		const auto physical_gorup_index = indexes[3];
		//const auto element_group_index	= value_set[4];

		auto& data = element_datas[i];

		data.figure = this->convert_to_figure(figure_type_index);
		data.figure_order = this->convert_to_figure_order(figure_type_index);
		data.element_type = physical_group_index_to_element_type.at(physical_gorup_index);

		constexpr auto num_type_indexes = 5;
		indexes.erase(indexes.begin(), indexes.begin() + num_type_indexes);

		const auto num_nodes = indexes.size();
		data.node_indexes.reserve(num_nodes);
		data.node_indexes.insert(data.node_indexes.end(), indexes.begin(), indexes.end());
	}

	LOG << std::left << std::setw(50) << "@ Extract " << num_elements << " Element Data " << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;

	return element_datas;
}

void Gmsh_Grid_File::write(const std::string_view output_file_path) const
{
	this->grid_text_.write(output_file_path);
}


ElementType Gmsh_Grid_File::convert_to_element_type(const Sentence& sentence) const
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
}

Figure Gmsh_Grid_File::convert_to_figure(const ushort figure_type_index) const
{
	switch (static_cast<Gmsh_Figure_Type>(figure_type_index))
	{
	case Gmsh_Figure_Type::POINT:		return Figure::point;
	case Gmsh_Figure_Type::LINE_P1:
	case Gmsh_Figure_Type::LINE_P2:
	case Gmsh_Figure_Type::LINE_P3:
	case Gmsh_Figure_Type::LINE_P4:
	case Gmsh_Figure_Type::LINE_P5:
	case Gmsh_Figure_Type::LINE_P6:		return Figure::line;
	case Gmsh_Figure_Type::TRIS_P1:
	case Gmsh_Figure_Type::TRIS_P2:
	case Gmsh_Figure_Type::TRIS_P3:
	case Gmsh_Figure_Type::TRIS_P4:
	case Gmsh_Figure_Type::TRIS_P5:		return Figure::triangle;
	case Gmsh_Figure_Type::QUAD_P1:
	case Gmsh_Figure_Type::QUAD_P2:
	case Gmsh_Figure_Type::QUAD_P3:
	case Gmsh_Figure_Type::QUAD_P4:
	case Gmsh_Figure_Type::QUAD_P5:
	case Gmsh_Figure_Type::QUAD_P6:		return Figure::quadrilateral;
	case Gmsh_Figure_Type::TETS_P1:
	case Gmsh_Figure_Type::TETS_P2:
	case Gmsh_Figure_Type::TETS_P3:
	case Gmsh_Figure_Type::TETS_P4:
	case Gmsh_Figure_Type::TETS_P5:		return Figure::tetrahedral;
	case Gmsh_Figure_Type::HEXA_P1:
	case Gmsh_Figure_Type::HEXA_P2:
	case Gmsh_Figure_Type::HEXA_P3:
	case Gmsh_Figure_Type::HEXA_P4:
	case Gmsh_Figure_Type::HEXA_P5:		return Figure::hexahedral;
	case Gmsh_Figure_Type::PRIS_P1:
	case Gmsh_Figure_Type::PRIS_P2:
	case Gmsh_Figure_Type::PRIS_P3:
	case Gmsh_Figure_Type::PRIS_P4:
	case Gmsh_Figure_Type::PRIS_P5:		return Figure::prism;
	case Gmsh_Figure_Type::PYRA_P1:
	case Gmsh_Figure_Type::PYRA_P2:
	case Gmsh_Figure_Type::PYRA_P3:
	case Gmsh_Figure_Type::PYRA_P4:		return Figure::pyramid;
	default:
		EXCEPTION("invalid element type index");
		return Figure::not_in_list;
	}
}

short Gmsh_Grid_File::convert_to_figure_order(const ushort figure_type_index) const
{
	switch (static_cast<Gmsh_Figure_Type>(figure_type_index))
	{
	case Gmsh_Figure_Type::POINT:		return 0;
	case Gmsh_Figure_Type::LINE_P1:
	case Gmsh_Figure_Type::TRIS_P1:
	case Gmsh_Figure_Type::QUAD_P1:
	case Gmsh_Figure_Type::TETS_P1:
	case Gmsh_Figure_Type::HEXA_P1:
	case Gmsh_Figure_Type::PRIS_P1:
	case Gmsh_Figure_Type::PYRA_P1:		return 1;
	case Gmsh_Figure_Type::LINE_P2:
	case Gmsh_Figure_Type::TRIS_P2:
	case Gmsh_Figure_Type::QUAD_P2:
	case Gmsh_Figure_Type::TETS_P2:
	case Gmsh_Figure_Type::HEXA_P2:
	case Gmsh_Figure_Type::PRIS_P2:
	case Gmsh_Figure_Type::PYRA_P2:		return 2;
	case Gmsh_Figure_Type::LINE_P3:
	case Gmsh_Figure_Type::TRIS_P3:
	case Gmsh_Figure_Type::QUAD_P3:
	case Gmsh_Figure_Type::TETS_P3:
	case Gmsh_Figure_Type::HEXA_P3:
	case Gmsh_Figure_Type::PRIS_P3:
	case Gmsh_Figure_Type::PYRA_P3:		return 3;
	case Gmsh_Figure_Type::LINE_P4:
	case Gmsh_Figure_Type::TRIS_P4:
	case Gmsh_Figure_Type::QUAD_P4:
	case Gmsh_Figure_Type::TETS_P4:
	case Gmsh_Figure_Type::HEXA_P4:
	case Gmsh_Figure_Type::PRIS_P4:
	case Gmsh_Figure_Type::PYRA_P4:		return 4;
	case Gmsh_Figure_Type::LINE_P5:
	case Gmsh_Figure_Type::TRIS_P5:
	case Gmsh_Figure_Type::QUAD_P5:
	case Gmsh_Figure_Type::TETS_P5:
	case Gmsh_Figure_Type::HEXA_P5:
	case Gmsh_Figure_Type::PRIS_P5:		return 5;
	case Gmsh_Figure_Type::LINE_P6:
	case Gmsh_Figure_Type::QUAD_P6:		return 6;
	default:
		EXCEPTION("invalid element type index");
		return -1;
	}
}

std::unordered_set<size_t> Gmsh_Grid_File::find_bdry_node_index_set(void) const
{
	const auto element_datas = this->extract_element_datas();
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

std::pair<size_t, size_t> Gmsh_Grid_File::find_start_end_line_index_pair(const Gmsh_Data_List target_data) const
{
	const auto keyword = this->gmsh_data_to_keyword(target_data);
	const auto line_index = this->grid_text_.find_line_index_has_keyword(keyword);
	const auto num_datas_line_index = line_index + 1;

	const auto& sentece = this->grid_text_[num_datas_line_index];
	const auto num_datas = sentece.to_value<size_t>();

	const auto start_line_index = line_index + 2;
	const auto end_line_index = start_line_index + num_datas;

	return { start_line_index, end_line_index };
}

Text Gmsh_Grid_File::extract_text_about(const Gmsh_Data_List target_data) const
{
	const auto [start_line_index, end_line_index] = this->find_start_end_line_index_pair(target_data);
	return this->grid_text_.extract(start_line_index, end_line_index);
}



std::string Gmsh_Grid_File::gmsh_data_to_keyword(const Gmsh_Data_List target_data) const
{
	switch (target_data)
	{
	case Gmsh_Data_List::node_data:
		return "$Nodes";
	case Gmsh_Data_List::element_data:
		return "$Elements";
	case Gmsh_Data_List::physical_name_data:
		return "$PhysicalNames";
	default:
		EXCEPTION("Gmsh does not have this data");
		return "";
	}
}

std::map<ushort, ElementType> Gmsh_Grid_File::make_physical_group_index_to_element_type(void) const
{
	const auto physical_name_text = this->extract_text_about(Gmsh_Data_List::physical_name_data);

	std::map<ushort, ElementType> physical_group_index_to_element_type;

	for (const auto& physical_name_sentence : physical_name_text)
	{
		const char delimiter = ' ';
		const auto parsed_sentences = physical_name_sentence.parse_by(delimiter);

		//const auto space_dimension	= parsed_sentence_set[0].to_value<ushort>();
		const auto physical_group_index = parsed_sentences[1].to_value<ushort>();
		const auto name					= parsed_sentences[2].get_remove("\"");
		const auto element_type			= this->convert_to_element_type(name);

		physical_group_index_to_element_type.emplace(physical_group_index, element_type);
	}

	return physical_group_index_to_element_type;
}
