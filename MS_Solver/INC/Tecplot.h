#pragma once
#include "Governing_Equation.h"
#include "Grid_Builder.h"


enum class Post_File_Type {
	Grid, Solution
};

enum class Zone_Type {
	FETriangle = 2,
	FETetraheron = 4
};

class Tecplot
{
private:
	using This_ = Tecplot;
	
public:
	static inline bool post_condition_ = false;

private:
	static inline std::string path_;
	static inline ushort post_order_ = 0;
	static inline std::string grid_variables_str_;
	static inline std::string solution_variables_str_;
	static inline Zone_Type zone_type_;
	static inline size_t num_element_ = 0;
	static inline size_t num_node_ = 0;
	static inline const double* time_ptr_ = nullptr;
	static inline std::vector<size_t> num_post_points_;
	static inline std::map<std::string, std::vector<double>> additioinal_data_name_to_values_;

	//For HOM 
	static inline std::vector<Dynamic_Matrix> set_of_basis_post_points_;

private:
	Tecplot(void) = delete;

public:
	static void set_path(const std::string& path) { This_::path_ = path; };

	template <typename Governing_Equation>
	static void initialize(const ushort post_order);

	static void syncronize_time(const double& current_time) { This_::time_ptr_ = &current_time; };

	template <ushort space_dimension>
	static void post_grid(const std::vector<Element<space_dimension>>& cell_elements);

public:
	template <typename T>
	static void record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);

	static void record_cell_indexes(void);

	template <typename T>
	static void conditionally_record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);

	static void conditionally_record_cell_indexes(void);

public: //for FVM
	template <ushort num_equation>	
	static void post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");

	template <ushort num_equation>
	static void conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");

public:	//for HOM
	template <ushort space_dimension, typename Reconstruction_Method>
	static void initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method);

	template <ushort num_equation, ushort num_basis>
	static void post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");

	template <ushort num_equation, ushort num_basis>
	static void conditionally_post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");

private:
	static void write_binary_header(const Post_File_Type file_type, const std::string_view post_file_path);
	static void write_binary_grid_post_file(const std::vector<std::vector<double>>& coordinates, const std::vector<std::vector<int>>& connectivities);
	static void write_binary_solution_post_file(const std::vector<std::vector<double>>& post_solution_binary_datas, const std::string& comment = "");

private:
	static constexpr bool is_scalar_equation(const ushort num_equation) { return num_equation == 1; };
	static void reset(void);
	static std::vector<int> convert_to_binary_data(const std::string& str);	

	template <ushort num_equation>
	static std::vector<std::vector<double>> convert_to_binary_data(const std::vector<Euclidean_Vector<num_equation>>& solutions);

	template <typename T>
	static std::vector<T> convert_cell_data_to_post_point_data(const std::vector<T>& cell_datas);
};


//template definition part
template <typename Governing_Equation>
void Tecplot::initialize(const ushort post_order) {
	This_::post_order_ = post_order;

	// Binary output
	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		This_::grid_variables_str_ = "X,Y";
		This_::solution_variables_str_ = "q";
		This_::zone_type_ = Zone_Type::FETriangle;
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		This_::grid_variables_str_ = "X,Y";
		This_::solution_variables_str_ = "rho,rhou,rhov,rhoE,u,v,p";
		This_::zone_type_ = Zone_Type::FETriangle;
	}
};

template <ushort space_dimension, typename Reconstruction_Method>
void Tecplot::initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method) {
	const auto& cell_elements = grid.elements.cell_elements;
	const auto num_cell = cell_elements.size();
	This_::set_of_basis_post_points_.reserve(num_cell);

	for (uint i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		auto post_nodes = geometry.post_nodes(This_::post_order_);
		This_::set_of_basis_post_points_.push_back(reconstruct_method.calculate_basis_nodes(i, post_nodes));
	}
}

template <ushort space_dimension>
void Tecplot::post_grid(const std::vector<Element<space_dimension>>& cell_elements) {
	//post processing grid data	
	const auto num_cell = cell_elements.size();
	This_::num_post_points_.resize(num_cell);

	size_t connectivity_start_index = 0;	//BINARY start with 0
	//size_t connectivity_start_index = 1;	//ASCII start with 1											

	std::vector<std::vector<double>> coordinates(space_dimension);
	std::vector<std::vector<int>> connectivities;

	for (uint i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		const auto post_nodes = geometry.post_nodes(This_::post_order_);
		for (const auto& node : post_nodes) {
			for (ushort j = 0; j < space_dimension; ++j) 
				coordinates[j].push_back(node[j]);
		}

		const auto post_connectivities = geometry.reference_geometry_.post_connectivities(This_::post_order_, connectivity_start_index);

		for (const auto& connectivity : post_connectivities) {
			const auto num_point = connectivity.size();			
			std::vector<int> temp(num_point);

			for (uint j = 0; j < num_point; ++j)
				temp[j] = static_cast<int>(connectivity[j]);

			connectivities.push_back(std::move(temp));
		}

		const auto num_post_node = post_nodes.size();
		connectivity_start_index += num_post_node;
		This_::num_node_ += num_post_node;
		This_::num_element_ += post_connectivities.size();
		This_::num_post_points_[i] = num_post_node;
	}

	This_::write_binary_grid_post_file(coordinates, connectivities);
}

template <typename T>
void Tecplot::record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas) {
	const auto num_cell = This_::num_post_points_.size();
	dynamic_require(cell_datas.size() == num_cell, "number of variable should be same with number of cell");

	std::vector<double> converted_cell_datas;
	converted_cell_datas.reserve(This_::num_node_);

	for (uint i = 0; i < num_cell; ++i)
		converted_cell_datas.push_back(static_cast<double>(cell_datas[i]));

	auto post_point_datas = This_::convert_cell_data_to_post_point_data(converted_cell_datas);

	if (This_::additioinal_data_name_to_values_.find(data_name) == This_::additioinal_data_name_to_values_.end())
		This_::additioinal_data_name_to_values_.emplace(data_name, std::move(post_point_datas));
	else
		This_::additioinal_data_name_to_values_.at(data_name) = std::move(post_point_datas);
}

template <typename T>
void Tecplot::conditionally_record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas) {
	if (This_::post_condition_)
		This_::record_cell_variables(data_name, cell_datas);
}

template <ushort num_equation>
void Tecplot::post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {	
	//FVM solutions post processing
	const auto post_point_solutions = This_::convert_cell_data_to_post_point_data(solutions);
	auto post_solution_binary_datas = This_::convert_to_binary_data(post_point_solutions);

	const auto num_solution_data = post_solution_binary_datas.size();
	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
	const auto num_total_data = num_solution_data + num_additional_data;
	post_solution_binary_datas.reserve(num_total_data);

	for (auto& [name, datas] : This_::additioinal_data_name_to_values_)
		post_solution_binary_datas.push_back(std::move(datas));

	This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
}

template <ushort num_equation>
void Tecplot::conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {
	if (This_::post_condition_)
		This_::post_solution(solutions, comment);
}

template <ushort num_equation, ushort num_basis>
void Tecplot::post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
	//HOM solution coefficients post processing
	std::vector<Euclidean_Vector<num_equation>> post_point_solutions;
	post_point_solutions.reserve(This_::num_node_);

	const auto num_solution_coefficient = solution_coefficients.size();
	for (uint i = 0; i < num_solution_coefficient; ++i) {
		const auto solution_post_points = solution_coefficients[i] * This_::set_of_basis_post_points_[i];
		for (ushort j = 0; j < This_::num_post_points_[i]; ++j)
			post_point_solutions.push_back(solution_post_points.column<num_equation>(j));
	}

	auto post_solution_binary_datas = This_::convert_to_binary_data(post_point_solutions);

	const auto num_solution_data = post_solution_binary_datas.size();
	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
	const auto num_total_data = num_solution_data + num_additional_data;
	post_solution_binary_datas.reserve(num_total_data);

	for (auto& [name, datas] : This_::additioinal_data_name_to_values_)
		post_solution_binary_datas.push_back(std::move(datas));

	This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
}

template <ushort num_equation, ushort num_basis>
void Tecplot::conditionally_post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
	if (This_::post_condition_)
		This_::post_solution(solution_coefficients, comment);
}


template <ushort num_equation>
std::vector<std::vector<double>> Tecplot::convert_to_binary_data(const std::vector<Euclidean_Vector<num_equation>>& post_point_solutions) {
	std::vector<std::vector<double>> convert_to_binary_data;

	if constexpr (This_::is_scalar_equation(num_equation)) {
		convert_to_binary_data.resize(num_equation);

		for (auto& post_point_values : convert_to_binary_data)
			post_point_values.reserve(This_::num_node_);

		for (ushort i = 0; i < num_equation; ++i)
			for (uint j = 0; j < This_::num_node_; ++j)
				convert_to_binary_data[i].push_back(post_point_solutions[j][i]);
	}
	else {
		convert_to_binary_data.resize(2 * num_equation - 1);

		for (auto& post_point_values : convert_to_binary_data)
			post_point_values.reserve(This_::num_node_);

		for (uint i = 0; i < This_::num_node_; ++i) {
			const auto& cvariable = post_point_solutions[i];
			const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);

			for (size_t j = 0; j < num_equation; ++j)
				convert_to_binary_data[j].push_back(cvariable[j]);

			for (size_t j = 0; j < num_equation - 1; ++j)
				convert_to_binary_data[j + 4].push_back(pvariable[j]);
		}
	}

	return convert_to_binary_data;
}

template<typename T>
std::vector<T> Tecplot::convert_cell_data_to_post_point_data(const std::vector<T>& cell_datas) {
	const auto num_cell = This_::num_post_points_.size();
	dynamic_require(cell_datas.size() == num_cell, "number of cell data should be same with num cell");

	std::vector<T> post_point_datas;
	post_point_datas.reserve(This_::num_node_);

	for (uint i = 0; i < num_cell; ++i) {
		for (ushort j = 0; j < num_post_points_[i]; ++j)
			post_point_datas.push_back(cell_datas[i]);
	}

	return post_point_datas;
}



//ASCII

//template <typename Governing_Equation>
//void Post_Solution_Data::initialize(const ushort post_order) {
//	This_::post_order_ = post_order;

//	//if constexpr (ms::is_SCL_2D<Governing_Equation>) {
//	//	This_::grid_variable_str_ = "Variables = X Y";
//	//	This_::solution_variable_str_ = "Variables = q";
//	//	This_::zone_type_str_ = "ZoneType = FETriangle";
//	//}
//	//else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
//	//	This_::grid_variable_str_ = "Variables = X Y";
//	//	This_::solution_variable_str_ = "Variables = rho rhou rhov rhoE u v p";
//	//	This_::zone_type_str_ = "ZoneType = FETriangle";
//	//}
//};

//Text Post_Solution_Data::header_text(const Post_File_Type file_type) {
//	static size_t strand_id = 0;
//
//	Text header;
//	header.reserve(10);
//	if (file_type == Post_File_Type::Grid) {
//		header << "Title = Grid";
//		header << "FileType = Grid";
//		header << This_::grid_variable_str_;
//		header << "Zone T = Grid";
//	}
//	else {
//		std::string solution_variable_str = This_::solution_variable_str_;
//		if (!This_::additioinal_data_name_to_values_.empty()) {
//			for (const auto& [info_name, info_values] : This_::additioinal_data_name_to_values_)
//				solution_variable_str += ", " + info_name;
//		}
//
//		header << "Title = Solution_at_" + ms::double_to_string(*time_ptr_);
//		header << "FileType = Solution";
//		header << solution_variable_str;
//		header << "Zone T = Solution_at_" + ms::double_to_string(*time_ptr_);
//
//		strand_id++;
//	}
//
//	header << This_::zone_type_str_;
//	header << "Nodes = " + std::to_string(num_node_);
//	header << "Elements = " + std::to_string(num_element_);
//	header << "DataPacking = Block";
//	header << "StrandID = " + std::to_string(strand_id);
//
//	if (file_type == Post_File_Type::Grid)
//		header << "SolutionTime = 0.0 \n\n";
//	else
//		header << "SolutionTime = " + ms::double_to_string(*time_ptr_) + "\n\n";
//
//	return header;
//}

//template <ushort space_dimension>
//void Post_Solution_Data::post_grid(const std::vector<Element<space_dimension>>& cell_elements) {

//	//const auto num_cell = cell_elements.size();
//	//This_::num_post_points_.resize(num_cell);
//
//	//ushort str_per_line = 1;
//	//size_t connectivity_start_index = 1;
//
//	//Text grid_post_data_text(space_dimension);
//	//for (uint i = 0; i < num_cell; ++i) {
//	//	const auto& geometry = cell_elements[i].geometry_;
//
//	//	const auto post_nodes = geometry.post_nodes(This_::post_order_);
//	//	for (const auto& node : post_nodes) {
//	//		for (ushort i = 0; i < space_dimension; ++i, ++str_per_line) {
//	//			grid_post_data_text[i] += ms::double_to_string(node.at(i)) + " ";
//	//			if (str_per_line == 10) {
//	//				grid_post_data_text[i] += "\n";
//	//				str_per_line = 1;
//	//			}
//	//		}
//	//	}
//
//	//	const auto connectivities = geometry.reference_geometry_.post_connectivities(This_::post_order_, connectivity_start_index);
//
//	//	std::string connectivity_str;
//	//	for (const auto& connectivity : connectivities) {
//	//		for (const auto index : connectivity)
//	//			connectivity_str += std::to_string(index) + " ";
//
//	//		grid_post_data_text << std::move(connectivity_str);
//	//	}
//
//	//	const auto num_post_node = post_nodes.size();
//	//	connectivity_start_index += num_post_node;
//	//	This_::num_node_ += num_post_node;
//	//	This_::num_element_ += connectivities.size();
//	//	This_::num_post_points_[i] = num_post_node;
//	//}
//
//	//auto grid_post_header_text = This_::header_text(Post_File_Type::Grid);
//
//	//const auto grid_file_path = This_::path_ + "grid.plt";
//	//grid_post_header_text.write(grid_file_path);
//	//grid_post_data_text.add_write(grid_file_path);
//}


//template <ushort num_equation>
//void Post_Solution_Data::write_solution_post_file(const std::vector<Euclidean_Vector<num_equation>>& post_point_solutions, const std::string& comment) {
//	static size_t count = 1;
//
//	std::string solution_file_path;
//	if (comment.empty())
//		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + ".plt";
//	else
//		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + "_" + comment + ".plt";
//
//	//solution post header text
//	//auto solution_post_header_text = This_::header_text(Post_File_Type::Solution);
//	//solution_post_header_text.write(solution_file_path);
//
//
//	//solution post data text
//	size_t str_per_line = 1;
//
//	Text solution_post_data_text;
//
//	if constexpr (This_::is_scalar_equation(num_equation)) {
//		solution_post_data_text.resize(num_equation);
//
//		for (size_t i = 0; i < This_::num_node_; ++i, ++str_per_line) {
//			const auto& solution = post_point_solutions[i];
//			solution_post_data_text[0] += ms::double_to_string(solution.at(0)) + " ";
//			if (str_per_line == 10) {
//				solution_post_data_text[0] += "\n";
//				str_per_line = 1;
//			}
//		}
//	}
//	else {
//		solution_post_data_text.resize(2 * num_equation);
//
//		for (size_t i = 0; i < This_::num_node_; ++i, ++str_per_line) {
//			const auto& cvariable = post_point_solutions[i];
//			const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);
//
//			//write conservative variable
//			for (size_t k = 0; k < num_equation; ++k)
//				solution_post_data_text[k] += ms::double_to_string(cvariable.at(k)) + " ";
//
//			//write primitive variable without a
//			for (size_t k = 0; k < num_equation - 1; ++k)
//				solution_post_data_text[k + 4] += ms::double_to_string(pvariable.at(k)) + " ";
//
//			if (str_per_line == 10) {
//				for (auto& sentence : solution_post_data_text)
//					sentence += "\n";
//
//				str_per_line = 1;
//			}
//		}
//	}
//
//	//additional data
//	str_per_line = 1;
//
//	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
//	Text additional_data_text;
//	additional_data_text.reserve(num_additional_data);
//
//	std::string data_str;
//	for (const auto& [name, values] : This_::additioinal_data_name_to_values_) {
//		for (size_t j = 0; j < This_::num_node_; ++j, ++str_per_line) {
//			data_str += ms::double_to_string(values[j]) + " ";
//
//			if (str_per_line == 10) {
//				data_str += "\n";
//
//				str_per_line = 1;
//			}
//		}
//
//		additional_data_text << std::move(data_str);
//	}
//
//	//merge & wirte
//	solution_post_data_text.merge(std::move(additional_data_text));
//	solution_post_data_text.add_write(solution_file_path);
//
//
//	This_::additioinal_data_name_to_values_.clear();
//	This_::is_time_to_post_ = false;
//	if (comment == "final") {
//		count = 1;
//		This_::reset();
//	}
//}