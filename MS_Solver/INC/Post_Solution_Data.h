#pragma once
#include "Governing_Equation.h"
#include "Grid_Builder.h"


enum class Post_File_Type {
	Grid, Solution
};


class Post_Solution_Data
{
private:
	using This_ = Post_Solution_Data;
	
public:
	static inline bool is_time_to_post_ = false;

private:
	static inline std::string path_;
	static inline ushort post_order_ = 0;
	static inline std::string grid_variable_str_;
	static inline std::string solution_variable_str_;
	static inline std::string zone_type_str_;
	static inline size_t num_element_ = 0;
	static inline size_t num_node_ = 0;
	static inline const double* time_ptr_ = nullptr;
	static inline std::vector<size_t> num_post_points_;
	static inline std::map<std::string, std::vector<double>> additioinal_data_name_to_values_;

	//For HOM 
	static inline std::vector<Dynamic_Matrix> set_of_basis_post_points_;

private:
	Post_Solution_Data(void) = delete;

public:
	static void set_path(const std::string& path) { This_::path_ = path; };

	template <typename Governing_Equation>
	static void initialize(const ushort post_order);

	static void syncronize_time(const double& current_time) { This_::time_ptr_ = &current_time; };

	template <ushort space_dimension>
	static void post_grid(const std::vector<Element<space_dimension>>& cell_elements);

	template <typename T>
	static void conditionally_record_cell_variables(const std::string& variable_name, const std::vector<T>& variables);

	static void conditionally_record_cell_indexes(void);

	template <typename T>
	static void record_cell_variables(const std::string& variable_name, const std::vector<T>& variables);

	static void record_cell_indexes(void);


public: //for FVM
	template <ushort num_equation>
	static void conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");

	template <ushort num_equation>	
	static void post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");



public:	//for HOM
	template <ushort space_dimension, typename Reconstruction_Method>
	static void initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method);

	template <ushort num_equation, ushort num_basis>
	static void conditionally_post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");

	template <ushort num_equation, ushort num_basis>
	static void post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");

private:
	static Text header_text(const Post_File_Type file_type);
	static constexpr bool is_scalar_equation(const ushort num_equation) { return num_equation == 1; };
	static void reset(void);

	template <ushort num_equation>
	static void write_solution_post_file(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");
};


//template definition part
template <typename Governing_Equation>
void Post_Solution_Data::initialize(const ushort post_order) {
	This_::post_order_ = post_order;

	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		This_::grid_variable_str_ = "Variables = X Y";
		This_::solution_variable_str_ = "Variables = q";
		This_::zone_type_str_ = "ZoneType = FETriangle";
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		This_::grid_variable_str_ = "Variables = X Y";
		This_::solution_variable_str_ = "Variables = rho rhou rhov rhoE u v p";
		This_::zone_type_str_ = "ZoneType = FETriangle";
	}
};

template <ushort space_dimension, typename Reconstruction_Method>
void Post_Solution_Data::initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method) {
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
void Post_Solution_Data::post_grid(const std::vector<Element<space_dimension>>& cell_elements) {
	const auto num_cell = cell_elements.size();
	This_::num_post_points_.resize(num_cell);

	ushort str_per_line = 1;
	size_t connectivity_start_index = 1;

	Text grid_post_data_text(space_dimension);
	for (uint i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		const auto post_nodes = geometry.post_nodes(This_::post_order_);
		for (const auto& node : post_nodes) {
			for (ushort i = 0; i < space_dimension; ++i, ++str_per_line) {
				grid_post_data_text[i] += ms::double_to_string(node.at(i)) + " ";
				if (str_per_line == 10) {
					grid_post_data_text[i] += "\n";
					str_per_line = 1;
				}
			}
		}

		const auto connectivities = geometry.reference_geometry_.post_connectivities(This_::post_order_, connectivity_start_index);

		std::string connectivity_str;
		for (const auto& connectivity : connectivities) {
			for (const auto index : connectivity)
				connectivity_str += std::to_string(index) + " ";

			grid_post_data_text << std::move(connectivity_str);
		}

		const auto num_post_node = post_nodes.size();
		connectivity_start_index += num_post_node;
		This_::num_node_ += num_post_node;
		This_::num_element_ += connectivities.size();
		This_::num_post_points_[i] = num_post_node;
	}

	auto grid_post_header_text = This_::header_text(Post_File_Type::Grid);

	const auto grid_file_path = This_::path_ + "grid.plt";
	grid_post_header_text.write(grid_file_path);
	grid_post_data_text.add_write(grid_file_path);
}


template <typename T>
void Post_Solution_Data::record_cell_variables(const std::string& variable_name, const std::vector<T>& variables) {
	const auto num_cell = This_::num_post_points_.size();
	dynamic_require(variables.size() == num_cell, "number of variable should be same with number of cell");

	std::vector<double> converted_variables;
	converted_variables.reserve(This_::num_node_);

	for (uint i = 0; i < num_cell; ++i) {
		const auto converted_variable = static_cast<double>(variables[i]);

		for (ushort j = 0; j < num_post_points_[i]; ++j)
			converted_variables.push_back(converted_variable);
	}

	This_::additioinal_data_name_to_values_.emplace(variable_name, std::move(converted_variables));
}


void Post_Solution_Data::record_cell_indexes(void) {
	const auto num_cell = num_post_points_.size();
	std::vector<uint> cell_index(num_cell);
	for (uint i = 0; i < num_cell; ++i)
		cell_index[i] = i;

	This_::record_cell_variables("cell_index", cell_index);
}

template <typename T>
void Post_Solution_Data::conditionally_record_cell_variables(const std::string& variable_name, const std::vector<T>& variables) {
	if (This_::is_time_to_post_)
		This_::record_cell_variables(variable_name, variables);
}

void Post_Solution_Data::conditionally_record_cell_indexes(void) {
	if (This_::is_time_to_post_)
		This_::record_cell_indexes();
}


template <ushort num_equation>
void Post_Solution_Data::conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {
	if (This_::is_time_to_post_)
		This_::post_solution(solutions, comment);
}

template <ushort num_equation>
void Post_Solution_Data::post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {	
	//FVM solutions post processing
	std::vector<Euclidean_Vector<num_equation>> post_point_solutions;
	post_point_solutions.reserve(This_::num_node_);

	const auto num_solution = solutions.size();

	for (uint i = 0; i < num_solution; ++i) {
		for (ushort j = 0; j < This_::num_post_points_[i]; ++j)
			post_point_solutions.push_back(solutions[i]);
	}

	This_::write_solution_post_file(post_point_solutions, comment);
}

template <ushort num_equation, ushort num_basis>
void Post_Solution_Data::conditionally_post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
	if (This_::is_time_to_post_)
		This_::post_solution(solution_coefficients, comment);
}

template <ushort num_equation, ushort num_basis>
void Post_Solution_Data::post_solution(const std::vector<Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
	//HOM solution coefficients post processing
	std::vector<Euclidean_Vector<num_equation>> post_point_solutions;
	post_point_solutions.reserve(This_::num_node_);

	const auto num_solution_coefficient = solution_coefficients.size();
	for (uint i = 0; i < num_solution_coefficient; ++i) {
		const auto solution_post_points = solution_coefficients[i] * This_::set_of_basis_post_points_[i];
		for (ushort j = 0; j < This_::num_post_points_[i]; ++j)
			post_point_solutions.push_back(solution_post_points.column<num_equation>(j));
	}

	This_::write_solution_post_file(post_point_solutions, comment);
}





Text Post_Solution_Data::header_text(const Post_File_Type file_type) {
	static size_t strand_id = 0;

	Text header;
	header.reserve(10);
	if (file_type == Post_File_Type::Grid) {
		header << "Title = Grid";
		header << "FileType = Grid";
		header << This_::grid_variable_str_;
		header << "Zone T = Grid";
	}
	else {
		std::string solution_variable_str = This_::solution_variable_str_;
		if (!This_::additioinal_data_name_to_values_.empty()) {
			for (const auto& [info_name, info_values] : This_::additioinal_data_name_to_values_)
				solution_variable_str += ", " + info_name;
		}

		header << "Title = Solution_at_" + ms::double_to_string(*time_ptr_);
		header << "FileType = Solution";
		header << solution_variable_str;
		header << "Zone T = Solution_at_" + ms::double_to_string(*time_ptr_);

		strand_id++;
	}

	header << This_::zone_type_str_;
	header << "Nodes = " + std::to_string(num_node_);
	header << "Elements = " + std::to_string(num_element_);
	header << "DataPacking = Block";
	header << "StrandID = " + std::to_string(strand_id);

	if (file_type == Post_File_Type::Grid)
		header << "SolutionTime = 0.0 \n\n";
	else
		header << "SolutionTime = " + ms::double_to_string(*time_ptr_) + "\n\n";

	return header;
}

void Post_Solution_Data::reset(void) {
	This_::is_time_to_post_ = false;
	This_::path_.clear();
	This_::post_order_ = 0;
	This_::grid_variable_str_.clear();
	This_::solution_variable_str_.clear();
	This_::zone_type_str_.clear();
	This_::num_element_ = 0;
	This_::num_node_ = 0;
	This_::time_ptr_ = nullptr;
	This_::num_post_points_.clear();
	This_::set_of_basis_post_points_.clear();
}



template <ushort num_equation>
void Post_Solution_Data::write_solution_post_file(const std::vector<Euclidean_Vector<num_equation>>& post_point_solutions, const std::string& comment) {
	static size_t count = 1;

	std::string solution_file_path;
	if (comment.empty())
		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + ".plt";
	else
		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + "_" + comment + ".plt";

	//solution post header text
	auto solution_post_header_text = This_::header_text(Post_File_Type::Solution);
	solution_post_header_text.write(solution_file_path);


	//solution post data text
	size_t str_per_line = 1;

	Text solution_post_data_text;

	if constexpr (This_::is_scalar_equation(num_equation)) {
		solution_post_data_text.resize(num_equation);

		for (size_t i = 0; i < This_::num_node_; ++i, ++str_per_line) {
			const auto& solution = post_point_solutions[i];
			solution_post_data_text[0] += ms::double_to_string(solution.at(0)) + " ";
			if (str_per_line == 10) {
				solution_post_data_text[0] += "\n";
				str_per_line = 1;
			}
		}
	}
	else {
		solution_post_data_text.resize(2 * num_equation);

		for (size_t i = 0; i < This_::num_node_; ++i, ++str_per_line) {
			const auto& cvariable = post_point_solutions[i];
			const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);

			//write conservative variable
			for (size_t k = 0; k < num_equation; ++k)
				solution_post_data_text[k] += ms::double_to_string(cvariable.at(k)) + " ";

			//write primitive variable
			for (size_t k = 0; k < num_equation; ++k)
				solution_post_data_text[k + 4] += ms::double_to_string(pvariable.at(k)) + " ";

			if (str_per_line == 10) {
				for (auto& sentence : solution_post_data_text)
					sentence += "\n";

				str_per_line = 1;
			}
		}
	}

	//additional data
	str_per_line = 1;
	
	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
	Text additional_data_text;
	additional_data_text.reserve(num_additional_data);

	std::string data_str;
	for (const auto& [name, values] : This_::additioinal_data_name_to_values_) {
		for (size_t j = 0; j < This_::num_node_; ++j, ++str_per_line) {
			data_str += ms::double_to_string(values[j]) + " ";

			if (str_per_line == 10) {
					data_str += "\n";
					 
				str_per_line = 1;
			}
		}

		additional_data_text << std::move(data_str);
	}

	//merge & wirte
	solution_post_data_text.merge(std::move(additional_data_text));
	solution_post_data_text.add_write(solution_file_path);


	This_::additioinal_data_name_to_values_.clear();
	This_::is_time_to_post_ = false;
	if (comment == "final") {
		count = 1;
		This_::reset();
	}
}