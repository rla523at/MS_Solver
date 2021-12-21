#pragma once
#include "Text.h"
#include "Post_Variables.h"
#include "Post_Variable_Convertor.h"

class Solution_Time_For_Header
{
public:
	virtual double compute(const double solution_time, const size_t strand_id) const
	{
		return solution_time;
	};
};

class Solution_Time_For_Debug : public Solution_Time_For_Header
{
public:
	double compute(const double solution_time, const size_t strand_id) const override
	{
		return 0.001 * strand_id;
	};
};

class Solution_Time_For_Header_Factory
{
public:
	static std::unique_ptr<Solution_Time_For_Header> make_unqiue(const Configuration& configuration)
	{
		const auto post_for_debug = configuration.get_post_for_debug();
		if (ms::compare_icase(post_for_debug, "Yes"))
		{
			return std::make_unique<Solution_Time_For_Debug>();
		}
		else
		{
			return std::make_unique<Solution_Time_For_Header>();
		}
	}
};

class Tecplot_File_Writer
{
public:
	Tecplot_File_Writer(std::unique_ptr<Solution_Time_For_Header>&& solution_time_for_header)
		: header_solution_time_(std::move(solution_time_for_header)) {};

public://Command
	void write_grid_file(const Post_Variables& post_variables, const std::string_view post_file_path);
	void write_solution_file(Post_Variables& post_variables, const std::string_view post_file_path);
	size_t strand_id(void) const;

protected:
	void set_common_header_variable(const Post_Variables& post_variables);
	virtual void set_grid_header_variable(const Post_Variables& post_variables) abstract;
	virtual void set_solution_header_variable(const Post_Variables& post_variables) abstract;
	virtual void write_header(const std::string_view post_file_path) abstract;
	virtual void write_grid_data(const Post_Variables& post_variables, const std::string_view post_file_path) const abstract;
	virtual void write_solution_data(const Post_Variables& post_variables, const std::string_view post_file_path) const abstract;
	
protected:
	//common header variable
	Zone_Type zone_type_ = Zone_Type::Noting;
	int num_post_points_ = 0;
	int num_post_elements_ = 0;
	double solution_time_ = 0.0;
	size_t strand_id_ = 0;

	std::unique_ptr<Solution_Time_For_Header> header_solution_time_;
};

class Tecplot_ASCII_File_Writer : public Tecplot_File_Writer
{
public:
	Tecplot_ASCII_File_Writer(std::unique_ptr<Solution_Time_For_Header>&& solution_time_for_header)
		: Tecplot_File_Writer(std::move(solution_time_for_header)) {};

private:
    void set_grid_header_variable(const Post_Variables& post_variables) override;
	void set_solution_header_variable(const Post_Variables& post_variables) override;
	void write_header(const std::string_view post_file_path) override;
    void write_grid_data(const Post_Variables& post_variables, const std::string_view post_file_path) const override;
	void write_solution_data(const Post_Variables& post_variables, const std::string_view post_file_path) const override;

	std::vector<std::vector<int>> convert_to_ASCII_connectivities(const std::vector<std::vector<int>>& connectivities) const;
	std::string make_post_variable_str(const std::vector<std::string>& post_variable_names) const;
	template <typename T>	void write_data(const std::vector<std::vector<T>>& set_of_post_datas, const std::string_view post_file_path) const  
	{
		ushort str_per_line = 0;
		const auto num_data = set_of_post_datas.size();

		Text ASCII_data_text;
		ASCII_data_text.add_empty_lines(num_data);

		for (int i = 0; i < num_data; ++i) 
		{
			const auto& post_datas = set_of_post_datas[i];

			for (const auto post_data : post_datas) 
			{
				ASCII_data_text[i].insert_with_space(post_data);
				if (++str_per_line == 10) 
				{
					ASCII_data_text[i] << "\n";
					str_per_line = 0;
				}
			}
			str_per_line = 0;
			ASCII_data_text[i] << "\n\n";
		}

		ASCII_data_text.add_write(post_file_path);
	}

private:
    //header variable for ASCII
    std::string title_;
    std::string file_type_str_;
    std::string post_variable_names_;
    std::string zone_title_;
    std::string variable_location_str_;
};

class Tecplot_Binary_File_Writer : public Tecplot_File_Writer
{
public:
	Tecplot_Binary_File_Writer(std::unique_ptr<Solution_Time_For_Header>&& solution_time_for_header)
		: Tecplot_File_Writer(std::move(solution_time_for_header)) {};

private:
	void set_grid_header_variable(const Post_Variables& post_variables) override;
	void set_solution_header_variable(const Post_Variables& post_variables) override;
	void write_header(const std::string_view post_file_path) override;
	void write_grid_data(const Post_Variables& post_variables, const std::string_view post_file_path) const override;
	void write_solution_data(const Post_Variables& post_variables, const std::string_view post_file_path) const override;

	std::vector<int> to_tecplot_binary_format(const std::string& str) const;
	std::vector<int> to_tecplot_binary_format(const std::vector<std::string>& strs) const;
	void write_data(const std::vector<std::vector<double>>& set_of_post_datas, const std::string_view post_file_path) const;

private:
	//header variable for Binary
	int file_type_ = -1;
	std::vector<int> title_tecplot_binary_format_;
	int num_variable_ = -1;
	std::vector<int> variable_names_tecplot_binary_format_;
	std::vector<int> zone_name_tecplot_binary_format_;
	int specify_variable_location_ = -1;
};

class Tecplot_File_Writer_Factory//static class
{
public:
	static std::unique_ptr<Tecplot_File_Writer> make_unique(const Configuration& configuration);

private:
	Tecplot_File_Writer_Factory(void) = delete;
};

namespace ms
{
	std::string to_string(const Zone_Type get_zone_type);
}











































//
//
//
//
//enum class Post_File_Type {
//	grid,
//	solution
//};
//
//
//enum class Post_File_Format {
//	ASCII,
//	binary
//};
//
//
//
//
//
//
//class Tecplot
//{
//private:
//	Tecplot(void) = delete;
//
//private:
//	using This_ = Tecplot;
//	
//public:
//	static inline bool post_condition_ = false;
//
//private:
//	static inline Post_File_Format file_format_ = Post_File_Format::binary;
//	static inline std::string path_;
//	static inline ushort post_order_ = 0;
//	static inline std::string grid_variables_str_;
//	static inline std::string solution_variables_str_;
//	static inline Zone_Type zone_type_;
//	static inline size_t num_element_ = 0;
//	static inline size_t num_node_ = 0;
//	static inline const double* time_ptr_ = nullptr;
//	static inline std::vector<size_t> num_post_points_;
//	static inline std::map<std::string, std::vector<double>> additioinal_data_name_to_values_;
//
//	//For HOM 
//	static inline std::vector<Matrix> set_of_basis_post_points_;
//
//public:
//	static void set_path(const std::string& path) { This_::path_ = path; };
//	static void syncronize_solution_time(const double& current_time) { This_::time_ptr_ = &current_time; };
//
//public:
//	template <typename T>
//	static void record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);
//
//	static void record_cell_indexes(void);
//
//	template <typename T>
//	static void conditionally_record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);
//
//	static void conditionally_record_cell_indexes(void);
//
//public:
//	template <typename Governing_Equation>
//	static void initialize(const ushort post_order, const Post_File_Format format);
//
//	template <ushort space_dimension>
//	static void post_grid(const Grid<space_dimension>& grid);
//
//public: //for FVM
//	template <ushort num_equation>	
//	static void post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");
//
//	template <ushort num_equation>
//	static void conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment = "");
//
//public:	//for HOM
//	template <ushort space_dimension, typename Reconstruction_Method>
//	static void initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method);
//
//	template <ushort num_equation, ushort num_basis>
//	static void post_solution(const std::vector<Static_Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");
//
//	template <ushort num_equation, ushort num_basis>
//	static void conditionally_post_solution(const std::vector<Static_Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment = "");
//
//private:
//	static void write_ASCII_header(const Post_File_Type file_type, const std::string_view post_file_path);
//	static void write_ASCII_grid_post_file(const std::vector<std::vector<double>>& post_coordinate_blocks, const std::vector<std::vector<int>>& connectivities);
//	static void write_ASCII_solution_post_file(const std::vector<std::vector<double>>& post_solution_datas, const std::string& comment = "");
//
//
//	static void write_binary_header(const Post_File_Type file_type, const std::string_view post_file_path);
//	static void write_binary_grid_post_file(const std::vector<std::vector<double>>& post_coordinate_blocks, const std::vector<std::vector<int>>& connectivities);
//	static void write_binary_solution_post_file(const std::vector<std::vector<double>>& post_solution_binary_datas, const std::string& comment = "");
//
//private:
//	static constexpr bool is_scalar_equation(const ushort num_equation) { return num_equation == 1; };
//	static void reset(void);
//	static std::vector<int> convert_to_binary_data(const std::string& str);	
//
//	template <ushort num_equation>
//	static std::vector<std::vector<double>> convert_to_binary_data(const std::vector<Euclidean_Vector<num_equation>>& solutions);
//
//	template <typename T>
//	static std::vector<T> convert_cell_data_to_post_point_data(const std::vector<T>& cell_datas);
//};
//
//namespace ms {
//	std::string to_string(const Zone_Type zone_type);
//}
//
////template definition part
//template <typename Governing_Equation>
//void Tecplot::initialize(const ushort post_order, const Post_File_Format format) {
//	This_::file_format_ = format;
//	This_::post_order_ = post_order;
//
//	if constexpr (Governing_Equation::space_dimension() == 2) {
//		if constexpr (ms::is_SCL<Governing_Equation>) {
//			This_::grid_variables_str_ = "X,Y";
//			This_::solution_variables_str_ = "q";
//			This_::zone_type_ = Zone_Type::FETriangle;
//		}
//		else if constexpr (ms::is_Euler<Governing_Equation>) {
//			This_::grid_variables_str_ = "X,Y";
//			This_::solution_variables_str_ = "rho,rhou,rhov,rhoE,u,v,p";
//			This_::zone_type_ = Zone_Type::FETriangle;
//		}		
//	}
//	else if constexpr (Governing_Equation::space_dimension() == 3) {
//		if constexpr (ms::is_SCL<Governing_Equation>) {
//			This_::grid_variables_str_ = "X,Y,Z";
//			This_::solution_variables_str_ = "q";
//			This_::zone_type_ = Zone_Type::FETetrahedron;
//		}
//		else if constexpr (ms::is_Euler<Governing_Equation>) {
//			This_::grid_variables_str_ = "X,Y,Z";
//			This_::solution_variables_str_ = "rho,rhou,rhov,rhow,rhoE,u,v,w,p";
//			This_::zone_type_ = Zone_Type::FETetrahedron;
//		}		
//	}
//};
//
//template <ushort space_dimension, typename Reconstruction_Method>
//void Tecplot::initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method) {
//	const auto set_of_post_nodes = grid.cell_set_of_post_nodes(This_::post_order_);
//	
//	const auto num_cell = set_of_post_nodes.size();
//	This_::set_of_basis_post_points_.reserve(num_cell);
//
//	for (uint i = 0; i < num_cell; ++i) 
//		This_::set_of_basis_post_points_.push_back(reconstruct_method.basis_nodes(i, set_of_post_nodes[i]));
//}
//
//template <ushort space_dimension>
//void Tecplot::post_grid(const Grid<space_dimension>& grid) {
//	const auto set_of_post_nodes = grid.cell_set_of_post_nodes(This_::post_order_);
//	
//	const auto num_cell = set_of_post_nodes.size();
//	This_::num_post_points_.reserve(num_cell);
//	for (const auto& post_nodes : set_of_post_nodes)
//		This_::num_post_points_.push_back(post_nodes.size());
//
//	const auto post_coordinate_blocks = grid.cell_post_coordinate_blocks(set_of_post_nodes);
//	const auto connectivities = grid.cell_set_of_connectivities(post_order_, set_of_post_nodes);
//
//	This_::num_node_ = post_coordinate_blocks.front().size();
//	This_::num_element_ = connectivities.size();
//
//	if (This_::file_format_ == Post_File_Format::binary)
//		This_::write_binary_grid_post_file(post_coordinate_blocks, connectivities);
//	else
//		This_::write_ASCII_grid_post_file(post_coordinate_blocks, connectivities);
//}
//
//template <typename T>
//void Tecplot::record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas) {
//	const auto num_cell = This_::num_post_points_.size();
//	dynamic_require(cell_datas.size() == num_cell, "number of variable should be same with number of cell");
//
//	std::vector<double> converted_cell_datas;
//	converted_cell_datas.reserve(This_::num_node_);
//
//	for (uint i = 0; i < num_cell; ++i)
//		converted_cell_datas.push_back(static_cast<double>(cell_datas[i]));
//
//	auto post_point_datas = This_::convert_cell_data_to_post_point_data(converted_cell_datas);
//
//	if (!This_::additioinal_data_name_to_values_.contains(data_name))
//		This_::additioinal_data_name_to_values_.emplace(data_name, std::move(post_point_datas));
//	else
//		This_::additioinal_data_name_to_values_.at(data_name) = std::move(post_point_datas);
//}
//
//template <typename T>
//void Tecplot::conditionally_record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas) {
//	if (This_::post_condition_)
//		This_::record_cell_variables(data_name, cell_datas);
//}
//
//template <ushort num_equation>
//void Tecplot::post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {	
//	//FVM solutions post processing
//	const auto post_point_solutions = This_::convert_cell_data_to_post_point_data(solutions);
//	auto post_solution_binary_datas = This_::convert_to_binary_data(post_point_solutions);
//
//	const auto num_solution_data = post_solution_binary_datas.size();
//	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
//	const auto num_total_data = num_solution_data + num_additional_data;
//	post_solution_binary_datas.reserve(num_total_data);
//
//	for (auto& [name, datas] : This_::additioinal_data_name_to_values_)
//		post_solution_binary_datas.push_back(std::move(datas));
//
//	if (This_::file_format_ == Post_File_Format::binary)
//		This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
//	else
//		This_::write_ASCII_solution_post_file(post_solution_binary_datas, comment);	
//}
//
//template <ushort num_equation>
//void Tecplot::conditionally_post_solution(const std::vector<Euclidean_Vector<num_equation>>& solutions, const std::string& comment) {
//	if (This_::post_condition_)
//		This_::post_solution(solutions, comment);
//}
//
//template <ushort num_equation, ushort num_basis>
//void Tecplot::post_solution(const std::vector<Static_Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
//	//HOM solution coefficients post processing
//	std::vector<Euclidean_Vector<num_equation>> post_point_solutions;
//	post_point_solutions.reserve(This_::num_node_);
//
//	const auto num_solution_coefficient = solution_coefficients.size();
//	for (uint i = 0; i < num_solution_coefficient; ++i) {
//		const auto solution_post_points = solution_coefficients[i] * This_::set_of_basis_post_points_[i];
//		for (ushort j = 0; j < This_::num_post_points_[i]; ++j)
//			post_point_solutions.push_back(solution_post_points.column<num_equation>(j));
//	}
//
//	auto post_solution_binary_datas = This_::convert_to_binary_data(post_point_solutions);
//
//	const auto num_solution_data = post_solution_binary_datas.size();
//	const auto num_additional_data = This_::additioinal_data_name_to_values_.size();
//	const auto num_total_data = num_solution_data + num_additional_data;
//	post_solution_binary_datas.reserve(num_total_data);
//
//	for (auto& [name, datas] : This_::additioinal_data_name_to_values_)
//		post_solution_binary_datas.push_back(std::move(datas));
//
//	if (This_::file_format_ == Post_File_Format::binary)
//		This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
//	else
//		This_::write_ASCII_solution_post_file(post_solution_binary_datas, comment);
//}
//
//template <ushort num_equation, ushort num_basis>
//void Tecplot::conditionally_post_solution(const std::vector<Static_Matrix<num_equation, num_basis>>& solution_coefficients, const std::string& comment) {
//	if (This_::post_condition_)
//		This_::post_solution(solution_coefficients, comment);
//}
//
//
//template <ushort num_equation>
//std::vector<std::vector<double>> Tecplot::convert_to_binary_data(const std::vector<Euclidean_Vector<num_equation>>& post_point_solutions) {
//	std::vector<std::vector<double>> convert_to_binary_data;
//
//	if constexpr (This_::is_scalar_equation(num_equation)) {
//		convert_to_binary_data.resize(num_equation);
//
//		for (auto& post_point_values : convert_to_binary_data)
//			post_point_values.reserve(This_::num_node_);
//
//		for (ushort i = 0; i < num_equation; ++i)
//			for (uint j = 0; j < This_::num_node_; ++j)
//				convert_to_binary_data[i].push_back(post_point_solutions[j][i]);
//	}
//	else {
//		//후에 추가적인 processig 부분이 다양하게 필요한 경우
//		//switch 역활을 하는 enum class 만들고 
//		//initialize 부분에서 Governing Equation에 따라서 switch 변수를 할당한뒤
//		//performance loss가 있지만 if 문으로 처리해야 Tecplot:: << 이 형식을 유지할 수 있다.
//
//		constexpr ushort space_dimension = num_equation - 2; //temporary code
//
//		convert_to_binary_data.resize(2 * num_equation - 1);
//
//		for (auto& post_point_values : convert_to_binary_data)
//			post_point_values.reserve(This_::num_node_);
//
//		for (uint i = 0; i < This_::num_node_; ++i) {
//			const auto& cvariable = post_point_solutions[i];
//			const auto pvariable = Euler<space_dimension>::conservative_to_primitive(cvariable);
//
//			for (size_t j = 0; j < num_equation; ++j)
//				convert_to_binary_data[j].push_back(cvariable[j]);
//
//			for (size_t j = 0; j < num_equation - 1; ++j)
//				convert_to_binary_data[num_equation + j].push_back(pvariable[j]);
//		}
//	}
//
//	return convert_to_binary_data;
//}
//
//template<typename T>
//std::vector<T> Tecplot::convert_cell_data_to_post_point_data(const std::vector<T>& cell_datas) {
//	const auto num_cell = This_::num_post_points_.size();
//	dynamic_require(cell_datas.size() == num_cell, "number of cell data should be same with num cell");
//
//	std::vector<T> post_point_datas;
//	post_point_datas.reserve(This_::num_node_);
//
//	for (uint i = 0; i < num_cell; ++i) {
//		for (ushort j = 0; j < num_post_points_[i]; ++j)
//			post_point_datas.push_back(cell_datas[i]);
//	}
//
//	return post_point_datas;
//}
