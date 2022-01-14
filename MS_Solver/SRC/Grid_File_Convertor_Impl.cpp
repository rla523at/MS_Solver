//#include "../INC/Grid_File_Convertor_Impl.h"
//
//Gmsh_Convertor::Gmsh_Convertor(const ushort space_dimension)
//{
//	this->space_dimension_ = space_dimension;
//}
//
//
//std::vector<Element> Gmsh_Convertor::convert_to_elements(const std::string_view grid_file_path) const
//{
//	Profiler::set_time_point();
//
//	Text grid_text;
//	grid_text.read(grid_file_path);
//
//	const auto node_text = this->gmsh_file_reader_.read_node_text(grid_text);
//	const auto element_text = this->gmsh_file_reader_.read_element_text(grid_text);
//	const auto physical_name_text = this->gmsh_file_reader_.read_physical_name_text(grid_text);
//
//	LOG << std::left << std::setw(50) << "@ Read Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;
//
//	Profiler::set_time_point();
//
//	const auto nodes = this->make_nodes(node_text);
//	const auto physical_group_index_to_element_type = this->make_physical_group_index_to_element_type(physical_name_text);
//	return this->make_elements(element_text, physical_group_index_to_element_type, nodes);
//	
//	//const auto elements = this->make_elements(element_text, physical_group_index_to_element_type, nodes);
//	//LOG << std::left << std::setw(50) << "@ Convert Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;
//	//return elements; //복사가 발생해서 error가 나는데 왜 어디서 복사가 발생하는거지 ?
//}
//
//std::vector<Euclidean_Vector> Gmsh_Convertor::make_nodes(const Text& node_text) const
//{
//	std::vector<Euclidean_Vector> node_datas;
//	node_datas.reserve(node_text.size());
//
//	for (const auto& node_sentence : node_text)
//	{
//		const char delimiter = ' ';
//		auto parsed_sentences = node_sentence.parse(delimiter);
//
//		std::vector<double> node_coords(this->space_dimension_);
//		for (ushort i = 0; i < this->space_dimension_; ++i)
//		{
//			node_coords[i] = parsed_sentences[i + 1].to_value<double>();// parsed_sentences[0] : node index
//		}
//
//		node_datas.push_back(std::move(node_coords));
//	}
//
//	return node_datas;
//}
//
//std::map<ushort, ElementType> Gmsh_Convertor::make_physical_group_index_to_element_type(const Text& physical_name_text) const
//{
//	std::map<ushort, ElementType> physical_group_index_to_element_type;
//
//	for (const auto& physical_name_sentence : physical_name_text)
//	{
//		const char delimiter = ' ';
//		const auto parsed_sentences = physical_name_sentence.parse(delimiter);
//
//		//const size_t dimension		= parsed_sentence_set[0].to_value<size_t>();
//		const auto physical_group_index = parsed_sentences[1].to_value<ushort>();
//		const auto name					= parsed_sentences[2].get_remove("\"");
//		const auto element_type			= ms::sentece_to_element_type(name);
//
//		physical_group_index_to_element_type.emplace(physical_group_index, element_type);
//	}
//
//	return physical_group_index_to_element_type;
//}
//
//
//std::vector<Element> Gmsh_Convertor::make_elements(const Text& element_text, const std::map<ushort, ElementType>& physical_group_index_to_element_type, const std::vector<Euclidean_Vector>& node_datas) const
//{
//	std::vector<Element> elements;
//	elements.reserve(element_text.size());
//
//	for (const auto& element_sentence : element_text)
//	{
//		const auto delimiter = ' ';
//		const auto parsed_sentences = element_sentence.parse(delimiter);
//
//		auto value_set = ms::sentences_to_values<uint>(parsed_sentences);
//
//		//const auto index					= value_set[0];
//		const auto figure_type_index		= value_set[1];
//		//const auto tag_index				= value_set[2];
//		const auto physical_gorup_index		= value_set[3];
//		//const auto element_group_index	= value_set[4];
//
//		//reference geometry
//		const auto figure = this->figure_type_index_to_figure(figure_type_index);
//		const auto figure_order = this->figure_type_index_to_figure_order(figure_type_index);
//		const auto& reference_geometry = Reference_Geometry_Container::get_shared_ptr(figure, figure_order);
//
//		//geometry
//		constexpr auto num_indexes = 5;
//		value_set.erase(value_set.begin(), value_set.begin() + num_indexes);
//
//		const auto num_nodes = value_set.size();
//		std::vector<uint> node_indexes(num_nodes);
//		for (ushort i = 0; i < num_nodes; ++i)
//		{
//			node_indexes[i] = value_set[i] - 1;		//Gmsh node index start with 1
//		}
//
//		auto nodes = ms::extract_by_index(node_datas, node_indexes);
//		Geometry geometry(reference_geometry, std::move(nodes));
//
//		//element
//		const auto type = physical_group_index_to_element_type.at(physical_gorup_index);
//		elements.push_back({ type, std::move(node_indexes), std::move(geometry) });
//	}
//	
//	LOG << std::left << std::setw(50) << "@ Convert Grid File" << " ----------- " << Profiler::get_time_duration() << "s\n\n" << Log::print_;
//
//	return elements;
//}
//
//std::unique_ptr<Grid_File_Convertor> Grid_File_Convertor_Factory::make_unique(const Configuration& configuration)
//{
//	const auto space_dimension = configuration.space_dimension();
//	const auto& grid_file_type = configuration.get_grid_file_type();
//
//	if (ms::contains_icase(grid_file_type, "gmsh"))
//	{
//		return std::make_unique<Gmsh_Convertor>(space_dimension);
//	}
//	else
//	{
//		EXCEPTION("grid file type in configuration file does not supported");
//		return nullptr;
//	}
//}
//
