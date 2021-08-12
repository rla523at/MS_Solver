#pragma once
#include "Governing_Equation.h"
#include "Spatial_Discrete_Method.h"
#include "Reconstruction_Method.h"
#include "Element.h"
#include "Text.h"


enum class Post_File_Type {
	Grid, Solution
};


template <typename Governing_Equation, ushort post_order>
class Post_Solution_Data_Base
{
private:
	static_require(ms::is_governing_equation<Governing_Equation>, "it should be governing equation");

	static constexpr size_t space_dimension_	= Governing_Equation::space_dimension();
	static constexpr size_t num_equation_		= Governing_Equation::num_equation();

protected:
	static inline std::string path_;
	static inline size_t num_element_ = 0;
	static inline size_t num_node_ = 0;
	static inline const double* time_ptr_ = nullptr;
	static inline std::vector<size_t> num_post_points_;

private:
	Post_Solution_Data_Base(void) = delete;

public:
	static void set_path(const std::string& path) { Post_Solution_Data_Base::path_ = path; };	
	static void syncronize_time(const double& current_time) { Post_Solution_Data_Base::time_ptr_ = &current_time; };

protected:
	static Text header_text(const Post_File_Type file_type);
};


template <typename Governing_Equation, ushort post_order>
class Post_FVM_Solution_Data : public Post_Solution_Data_Base<Governing_Equation, post_order>
{
private:
	static_require(ms::is_governing_equation<Governing_Equation>, "it should be governing equation");

	static constexpr size_t space_dimension_	= Governing_Equation::space_dimension();
	static constexpr size_t num_equation_		= Governing_Equation::num_equation();

	using This_ = Post_FVM_Solution_Data<Governing_Equation, post_order>;

private:
	Post_FVM_Solution_Data(void) = delete;

public:
	static void post_grid(const std::vector<Element<space_dimension_>>& cell_elements);
	static void post_solution(const std::vector<Euclidean_Vector<num_equation_>>& solutions, const std::string& comment = "");
};


template <typename Governing_Equation, typename Reconstruction_Method, ushort post_order>
class Post_HOM_Solution_Data : public Post_Solution_Data_Base<Governing_Equation, post_order>
{
private:
	static_require(ms::is_governing_equation<Governing_Equation>, "it should be governing equation");

	static constexpr ushort space_dimension_	= Governing_Equation::space_dimension();
	static constexpr ushort num_equation_		= Governing_Equation::num_equation();
	static constexpr ushort num_basis_			= Reconstruction_Method::num_basis();

	using This_ = Post_HOM_Solution_Data<Governing_Equation, Reconstruction_Method, post_order>;

private:
	static inline std::vector<Dynamic_Matrix> set_of_basis_post_points_;

public:
	static void post_grid(const std::vector<Element<space_dimension_>>& cell_elements);
	static void post_solution(const std::vector<Matrix<num_equation_,num_basis_>>& solution_coefficients, const std::string& comment = "");
};


template <typename Governing_Equation, typename Spatial_Discrete_Method, typename Reconstruction_Method, ushort post_order>
class Post_Solution_Data {};


template <typename Governing_Equation, typename Reconstruction_Method, ushort post_order>
class Post_Solution_Data<Governing_Equation, FVM, Reconstruction_Method, post_order> : public Post_FVM_Solution_Data<Governing_Equation, post_order>
{};

template <typename Governing_Equation, typename Reconstruction_Method, ushort post_order>
class Post_Solution_Data<Governing_Equation, HOM, Reconstruction_Method, post_order> : public Post_HOM_Solution_Data<Governing_Equation, Reconstruction_Method, post_order>
{};


//template definition part
template <typename Governing_Equation, ushort post_order>
Text Post_Solution_Data_Base<Governing_Equation, post_order>::header_text(const Post_File_Type file_type) {
	static size_t strand_id = 0;

	std::string grid_variable_str;
	std::string solution_variable_str;
	std::string zone_type_str;


	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		grid_variable_str		= "Variables = X Y";
		solution_variable_str	= "Variables = q";
		zone_type_str			= "ZoneType = FETriangle";
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		grid_variable_str		= "Variables = X Y";
		solution_variable_str	= "Variables = rho rhou rhov rhoE u v p a";
		zone_type_str			= "ZoneType = FETriangle";
	}


	Text header;
	header.reserve(10);
	if (file_type == Post_File_Type::Grid) {
		header << "Title = Grid";
		header << "FileType = Grid";
		header << grid_variable_str;
		header << "Zone T = Grid";
	}
	else {
		header << "Title = Solution_at_" + ms::double_to_string(*time_ptr_);
		header << "FileType = Solution";
		header << solution_variable_str;
		strand_id++;
		header << "Zone T = Solution_at_" + ms::double_to_string(*time_ptr_);
	}

	header << zone_type_str;
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


template <typename Governing_Equation, ushort post_order>
void Post_FVM_Solution_Data<Governing_Equation, post_order>::post_grid(const std::vector<Element<space_dimension_>>& cell_elements) {
	const auto num_cell = cell_elements.size();
	This_::num_post_points_.resize(num_cell);

	ushort str_per_line = 1;
	size_t connectivity_start_index = 1;

	Text grid_post_data_text(space_dimension_);
	for (uint i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		const auto post_nodes = geometry.post_nodes(post_order);
		for (const auto& node : post_nodes) {
			for (ushort i = 0; i < space_dimension_; ++i, ++str_per_line) {
				grid_post_data_text[i] += ms::double_to_string(node.at(i)) + " ";
				if (str_per_line == 10) {
					grid_post_data_text[i] += "\n";
					str_per_line = 1;
				}
			}
		}

		const auto connectivities = geometry.reference_geometry_.post_connectivities(post_order, connectivity_start_index);

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


template <typename Governing_Equation, ushort post_order>
void Post_FVM_Solution_Data<Governing_Equation, post_order>::post_solution(const std::vector<Euclidean_Vector<num_equation_>>& solutions, const std::string& comment) {
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

	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		solution_post_data_text.resize(num_equation_);

		const auto num_solution = solutions.size();
		for (size_t i = 0; i < num_solution; ++i) {
			const auto& solution = solutions[i];
			for (size_t j = 0; j < This_::num_post_points_[i]; ++j, ++str_per_line) {
				solution_post_data_text[0] += ms::double_to_string(solution.at(0)) + " ";
				if (str_per_line == 10) {
					solution_post_data_text[0] += "\n";
					str_per_line = 1;
				}
			}
		}
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		solution_post_data_text.resize(2 * num_equation_);

		const auto num_solution = solutions.size();
		for (size_t i = 0; i < num_solution; ++i) {
			const auto& cvariable = solutions[i];
			const auto pvariable = Euler_2D::conservative_to_primitive(cvariable);

			for (size_t j = 0; j < This_::num_post_points_[i]; ++j, ++str_per_line) {
				//write conservative variable
				for (size_t k = 0; k < This_::num_equation_; ++k)
					solution_post_data_text[k] += ms::double_to_string(cvariable.at(k)) + " ";

				//write primitive variable
				for (size_t k = 0; k < This_::num_equation_; ++k)
					solution_post_data_text[k + 4] += ms::double_to_string(pvariable.at(k)) + " ";

				if (str_per_line == 10) {
					for (auto& sentence : solution_post_data_text)
						sentence += "\n";

					str_per_line = 1;
				}
			}
		}
	}

	solution_post_data_text.add_write(solution_file_path);
}


template <typename Governing_Equation, typename Reconstruction_Method, ushort post_order>
void Post_HOM_Solution_Data<Governing_Equation, Reconstruction_Method, post_order>::post_grid(const std::vector<Element<space_dimension_>>& cell_elements) {
	const auto num_cell = cell_elements.size();

	std::vector<std::vector<Euclidean_Vector<space_dimension_>>> set_of_post_nodes;
	set_of_post_nodes.reserve(num_cell);
	
	ushort str_per_line = 1;
	size_t connectivity_start_index = 1;

	Text grid_post_data_text(space_dimension_);
	for (uint i = 0; i < num_cell; ++i) {
		const auto& geometry = cell_elements[i].geometry_;

		auto post_nodes = geometry.post_nodes(post_order);

		for (const auto& node : post_nodes) {
			for (ushort i = 0; i < space_dimension_; ++i, ++str_per_line) {
				grid_post_data_text[i] += ms::double_to_string(node.at(i)) + " ";
				if (str_per_line == 10) {
					grid_post_data_text[i] += "\n";
					str_per_line = 1;
				}
			}
		}
		
		const auto connectivities = geometry.reference_geometry_.post_connectivities(post_order, connectivity_start_index);

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

		set_of_post_nodes.push_back(std::move(post_nodes));
	}	

	This_::set_of_basis_post_points_ = Reconstruction_Method::calculate_set_of_basis_nodes(set_of_post_nodes);

	auto grid_post_header_text = This_::header_text(Post_File_Type::Grid);

	const auto grid_file_path = This_::path_ + "grid.plt";
	grid_post_header_text.write(grid_file_path);
	grid_post_data_text.add_write(grid_file_path);
}


template <typename Governing_Equation, typename Reconstruction_Method, ushort post_order>
void Post_HOM_Solution_Data<Governing_Equation, Reconstruction_Method, post_order>::post_solution(const std::vector<Matrix<num_equation_, num_basis_>>& solution_coefficients, const std::string& comment) {
	static ushort count = 1;

	std::string solution_file_path;
	if (comment.empty())
		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + ".plt";
	else
		solution_file_path = This_::path_ + "solution_" + std::to_string(count++) + "_" + comment + ".plt";

	//solution post header text
	auto solution_post_header_text = This_::header_text(Post_File_Type::Solution);
	solution_post_header_text.write(solution_file_path);

	//solution post data text
	ushort str_per_line = 1;

	Text solution_post_data_text;
	const auto num_solution = solution_coefficients.size();

	solution_post_data_text.resize(This_::num_equation_);

	if constexpr (ms::is_SCL_2D<Governing_Equation>) {
		for (uint i = 0; i < num_solution; ++i) {
			const auto post_point_solutions = solution_coefficients[i] * This_::set_of_basis_post_points_[i];

			for (ushort j = 0; j < This_::num_post_points_[i]; ++j, ++str_per_line) {
				const auto solution = post_point_solutions.column(j);
				solution_post_data_text[0] += ms::double_to_string(solution[0]) + " ";
				if (str_per_line == 10) {
					solution_post_data_text[0] += "\n";
					str_per_line = 1;
				}
			}
		}
	}
	else if constexpr (std::is_same_v<Governing_Equation, Euler_2D>) {
		solution_post_data_text.resize(2 * This_::num_equation_);

		for (size_t i = 0; i < num_solution; ++i) {
			const auto post_point_solutions = solution_coefficients[i] * This_::set_of_basis_post_points_[i];

			for (ushort j = 0; j < This_::num_post_points_[i]; ++j, ++str_per_line) {
				const auto cvariable = post_point_solutions.column<This_::num_equation_>(j);
				const auto pvariable = Governing_Equation::conservative_to_primitive(cvariable);

				//write conservative variable
				for (size_t k = 0; k < This_::num_equation_; ++k) 
					solution_post_data_text[k] += ms::double_to_string(cvariable.at(k)) + " ";

				//write primitive variable
				for (size_t k = 0; k < This_::num_equation_; ++k) 
					solution_post_data_text[k + 4] += ms::double_to_string(pvariable.at(k)) + " ";

				if (str_per_line == 10) {					
					for (auto& sentence : solution_post_data_text)
						sentence += "\n";

					str_per_line = 1;
				}
			}
		}
	}

	solution_post_data_text.add_write(solution_file_path);
}