#pragma once
#include "Governing_Equation.h"
#include "Grid.h"

namespace ms {
	template <typename T, typename... Args>
	std::vector<T>& merge(std::vector<T>& vec1, std::vector<T>&& vec2, Args&&... args) {
		static_require((... && std::is_same_v<Args, std::vector<T>>), "every arguments should be vector of same type");

		vec1.reserve(vec1.size() + vec2.size());
		vec1.insert(vec1.end(), std::make_move_iterator(vec2.begin()), std::make_move_iterator(vec2.end()));

		if constexpr (sizeof...(Args) == 0)
			return vec1;
		else
			return ms::merge(vec1, std::move(args)...);
	}
}

enum class Zone_Type {
	FETriangle = 2,
	FETetrahedron = 4
};

enum class Tecplot_Header_Writer_Mode {
	grid,
	solution
};

class Post_Variables
{
public:
	void record_grid(const Grid<2>& grid) const;
	void record_variable(const std::string_view name, const std::vector<double>& values);

	Zone_Type zone_type(void) const;
	size_t num_post_node(void) const;
	size_t num_element(void) const;
	ushort num_grid_variable(void) const;
	ushort num_solution_variable(void) const;

	std::string grid_variable_str(void) const;
	std::string solution_variable_str(void) const;
	std::string variable_location_str(void) const;

	std::vector<std::vector<double>> calculate_set_of_grid_datas(void) const;
	std::vector<std::vector<double>> calculate_set_of_ASCII_grid_datas(void) const {
		//auto set_of_grid_datas = post_variables.calculate_set_of_grid_datas();
		//const auto num_grid_variable = post_variables.num_grid_variable();
		//
		//for (size_t i = num_grid_variable; i < set_of_grid_datas.size(); ++i) {
		//	auto& grid_datas = set_of_grid_datas[i];
		//	for (auto& data : grid_datas)
		//		data += 1; //ASCII connectivity start with 1
		//}
	};
	std::vector<std::vector<double>> calculate_set_of_solution_datas(void) const;
};

class TecPlot_Header_Writer_Base
{
public:
	void syncronize_time(const double& solution_time) {
		this->solution_time_ptr_ = &solution_time;
	}
	void write_grid_header(const Post_Variables& post_variables, const std::string_view post_file_path) {
		this->grid_mode(post_variables);
		this->write_header(post_file_path);
	}
	void write_solution_header(const Post_Variables& post_variables, const std::string_view post_file_path) {
		this->solution_mode(post_variables);
		this->write_header(post_file_path);
	}

protected:
	virtual void grid_mode(const Post_Variables& post_variables) abstract;
	virtual void solution_mode(const Post_Variables& post_variables) abstract;
	virtual void write_header(const std::string_view post_file_path) const abstract;

protected:
	const double* solution_time_ptr_ = nullptr;
	Tecplot_Header_Writer_Mode mode_;
	Zone_Type zone_type_;
	int num_post_nodes_ = 0;
	int num_elements_ = 0;
};

class Tecplot_ASCII_Header_Writer : public TecPlot_Header_Writer_Base
{
private:
	void grid_mode(const Post_Variables& post_variables) override {
		this->zone_type_ = post_variables.zone_type();
		this->num_post_nodes_ = static_cast<int>(post_variables.num_post_node());
		this->num_elements_ = static_cast<int>(post_variables.num_element());
		this->title_ = "Grid";
		this->file_type_str_ = "Grid";
		this->variable_names_ = post_variables.grid_variable_str();
		this->zone_title_ = "Grid";
		this->variable_location_str_ = "()";

		this->mode_ = Tecplot_Header_Writer_Mode::grid;
	}
	void solution_mode(const Post_Variables& post_variables) override {
		this->zone_type_ = post_variables.zone_type();
		this->num_post_nodes_ = static_cast<int>(post_variables.num_post_node());
		this->num_elements_ = static_cast<int>(post_variables.num_element());
		this->title_ = "Solution_at_" + ms::double_to_string(*this->solution_time_ptr_);
		this->file_type_str_ = "Solution";
		this->variable_names_ = post_variables.solution_variable_str();
		this->zone_title_ = "Solution_at_" + ms::double_to_string(*this->solution_time_ptr_);
		this->variable_location_str_ = post_variables.variable_location_str();

		this->mode_ = Tecplot_Header_Writer_Mode::solution;
	}
	void write_header(const std::string_view post_file_path) const override {
		//nullptr 예외 발생 코드
		
		static size_t strand_id = 0;

		Text header;
		header.reserve(11);
		header << "Title = " + this->title_;
		header << "FileType = " + this->file_type_str_;
		header << "Variables = " + this->variable_names_;
		header << "Zone T = " + this->zone_title_;
		header << "ZoneType =  " + ms::to_string(this->zone_type_);
		header << "Nodes = " + std::to_string(this->num_post_nodes_);
		header << "Elements = " + std::to_string(this->num_elements_);
		header << "DataPacking = Block";
		header << "StrandID = " + std::to_string(strand_id);		
		header << "SolutionTime = " + std::to_string(*this->solution_time_ptr_);
		header << "VarLocation = " + variable_location_str_;

		header.write(post_file_path);
	}

protected:
	std::string title_;
	std::string file_type_str_;
	std::string variable_names_;
	std::string zone_title_;
	std::string variable_location_str_;
};

class Tecplot_Binary_Header_Writer : public TecPlot_Header_Writer_Base
{
public:
	void grid_mode(const Post_Variables& post_variables) override {
		this->zone_type_ = post_variables.zone_type();
		this->num_post_nodes_ = static_cast<int>(post_variables.num_post_node());
		this->num_elements_ = static_cast<int>(post_variables.num_element());
		this->file_type_ = 1;
		this->title_tecplot_binary_format_ = this->to_tecplot_binary_format("Grid");
		this->num_variable_ = static_cast<int>(post_variables.num_grid_variable());
		this->variable_names_tecplot_binary_format_ = this->to_tecplot_binary_format(ms::parse(post_variables.grid_variable_str(), ','));
		this->zone_name_tecplot_binary_format_ = this->to_tecplot_binary_format("Grid");
		this->specify_variable_location_ = 0;

		this->mode_ = Tecplot_Header_Writer_Mode::grid;
	}
	void solution_mode(const Post_Variables& post_variables) override {
		this->zone_type_ = post_variables.zone_type();
		this->num_post_nodes_ = static_cast<int>(post_variables.num_post_node());
		this->num_elements_ = static_cast<int>(post_variables.num_element());
		this->file_type_ = 2;
		this->title_tecplot_binary_format_ = this->to_tecplot_binary_format("Solution_at_" + std::to_string(*this->solution_time_ptr_));
		this->num_variable_ = static_cast<int>(post_variables.num_solution_variable());
		this->variable_names_tecplot_binary_format_ = this->to_tecplot_binary_format(ms::parse(post_variables.solution_variable_str(), ','));
		this->zone_name_tecplot_binary_format_ = this->to_tecplot_binary_format("Solution_at_" + std::to_string(*this->solution_time_ptr_));

		if (post_variables.variable_location_str() == "()")
			this->specify_variable_location_ = 0;
		else
			this->specify_variable_location_ = 1;

		this->mode_ = Tecplot_Header_Writer_Mode::solution;
	}
	void write_header(const std::string_view post_file_path) const override {
		static int strand_id = 0;

		Binary_Writer writer(post_file_path);

		//I. HEADER SECTION

		//i. version number
		writer << "#!TDV112";												//default	

		//ii
		writer << 1;														//byte order - default

		//iii. Title and variable names

		writer << this->file_type_;
		writer << this->title_tecplot_binary_format_;
		writer << this->num_variable_;
		writer << this->variable_names_tecplot_binary_format_;

		//iiii. zones
		writer << 299.0f;													//zone marker
		writer << this->zone_name_tecplot_binary_format_;					//zone name
		writer << -1;														//parent zone - default
		writer << strand_id;												//strand id
		writer << *this->solution_time_ptr_;								//solution_time
		writer << -1;														//not used - default
		writer << static_cast<int>(this->zone_type_);						//zone type
		if (this->specify_variable_location_ == 0)							//specify var location
			writer << 0;
		else
			writer << 1 << 1;
		writer << 0;														//one to one face neighbor - default
		writer << 0;														//user define face neighbor connection -default
		writer << this->num_post_nodes_;									//num points
		writer << this->num_elements_;										//num elements
		writer << 0 << 0 << 0;												//i,j,k cell dim - default
		writer << 0;														//auxilarily name index pair - default
		writer << 357.0f;													//EOH_marker
	}

private:
	std::vector<int> to_tecplot_binary_format(const std::string& str) const {
		std::vector<int> tecplot_binary_format;
		tecplot_binary_format.insert(tecplot_binary_format.end(), str.begin(), str.end());
		tecplot_binary_format.push_back(0); // null
		return tecplot_binary_format;
	}
	std::vector<int> to_tecplot_binary_format(const std::vector<std::string>& strs) const {
		std::vector<int> tecplot_binary_format;

		for (const auto& str : strs) {
			auto str_binary_format = this->to_tecplot_binary_format(str);
			ms::merge(tecplot_binary_format, std::move(str_binary_format));
		}

		return tecplot_binary_format;
	}

protected:
	int file_type_;
	std::vector<int> title_tecplot_binary_format_;
	int num_variable_;
	std::vector<int> variable_names_tecplot_binary_format_;
	std::vector<int> zone_name_tecplot_binary_format_;
	int specify_variable_location_;
};

class Tecplot_ASCII_Data_Writer
{
public:
	void write_grid_data(const Post_Variables& post_variables, const std::string_view post_file_path) const {
		this->write_data(post_variables.calculate_set_of_ASCII_grid_datas(), post_file_path);
	}
	void write_solution_data(const Post_Variables& post_variables, const std::string_view post_file_path) const {
		this->write_data(post_variables.calculate_set_of_solution_datas(), post_file_path);
	}

private:
	void write_data(const std::vector<std::vector<double>>& set_of_post_datas, const std::string_view post_file_path) const {
		ushort str_per_line = 0;
		const auto num_data = set_of_post_datas.size();

		Text ASCII_data_text(num_data);

		for (ushort i = 0; i < num_data; ++i) {
			const auto& post_datas = set_of_post_datas[i];
			for (const auto post_data : post_datas) {
				ASCII_data_text[i] += ms::double_to_string(post_data) + " ";
				if (++str_per_line == 10) {
					ASCII_data_text[i] += "\n";
					str_per_line = 0;
				}
				ASCII_data_text[i] += "\n\n\n";
			}
		}

		ASCII_data_text.add_write(post_file_path);
	};
};

class Tecplot_Binary_Data_Writer
{
public:
	void write_grid_data(const Post_Variables& post_variables, const std::string_view post_file_path) const {
		this->write_data(post_variables.calculate_set_of_grid_datas(), post_variables.num_grid_variable(), post_file_path);
	}
	void write_solution_data(const Post_Variables& post_variables, const std::string_view post_file_path) const {
		this->write_data(post_variables.calculate_set_of_solution_datas(), post_variables.num_solution_variable(), post_file_path);
	}

private:
	void write_data(const std::vector<std::vector<double>>& set_of_post_datas, const size_t num_post_variable, const std::string_view post_file_path) const {
		//const auto set_of_post_datas = post_data_builder_->calculate_set_of_post_datas(post_variables);	//data = post variable + connectivity
		//const auto num_post_variable = post_data_builder_->calculate_num_variable(post_variables);

		//II. DATA SECTION		
		Binary_Writer binary_data_file(post_file_path, std::ios::app);

		//zone
		binary_data_file << 299.0f;								//zone marker

		for (ushort i = 0; i < num_post_variable; ++i)
			binary_data_file << 2;								//variable data format, double = 2

		binary_data_file << 0 << 0 << -1;						//has passive variable, has variable sharing, zone number to share connectivity - default

		for (ushort i = 0; i < num_post_variable; ++i) {
			const auto min_value = *std::min_element(set_of_post_datas[i].begin(), set_of_post_datas[i].end());
			const auto max_value = *std::max_element(set_of_post_datas[i].begin(), set_of_post_datas[i].end());
			binary_data_file << min_value << max_value;			//min,max value of each variable
		}

		for (ushort i = 0; i < set_of_post_datas.size(); ++i)
			binary_data_file << set_of_post_datas[i];			//values of datas
	};
};

class Tecplot_File_Writer abstract
{
public:
	virtual void write_grid_file(const Post_Variables& post_variables, const std::string_view post_file_path) abstract;
	virtual void write_solution_file(const Post_Variables& post_variables, const std::string_view post_file_path) abstract;
};

class Tecplot_ASCII_File_Writer : public Tecplot_File_Writer
{
public:
	void write_grid_file(const Post_Variables& post_variables, const std::string_view post_file_path) override {
		tecplot_header_writer_.write_grid_header(post_variables, post_file_path);
		tecplot_data_writer_.write_grid_data(post_variables, post_file_path);
	}
	void write_solution_file(const Post_Variables& post_variables, const std::string_view post_file_path) override {
		tecplot_header_writer_.write_solution_header(post_variables, post_file_path);
		tecplot_data_writer_.write_solution_data(post_variables, post_file_path);
	}

private:
	Tecplot_ASCII_Header_Writer tecplot_header_writer_;
	Tecplot_ASCII_Data_Writer tecplot_data_writer_;
};

class Tecplot_Binary_File_Writer : public Tecplot_File_Writer
{
public:
	void write_grid_file(const Post_Variables& post_variables, const std::string_view post_file_path) override {
		tecplot_header_writer_.write_grid_header(post_variables, post_file_path);
		tecplot_data_writer_.write_grid_data(post_variables, post_file_path);
	}
	void write_solution_file(const Post_Variables& post_variables, const std::string_view post_file_path) override {
		tecplot_header_writer_.write_solution_header(post_variables, post_file_path);
		tecplot_data_writer_.write_solution_data(post_variables, post_file_path);
	}

private:
	Tecplot_Binary_Header_Writer tecplot_header_writer_;
	Tecplot_Binary_Data_Writer tecplot_data_writer_;
};


class Discretized_Solution 
{
public:
	std::vector<std::vector<double>> calculate_post_point_solutions_by_variables(void) const;
	
	const std::vector<std::string_view>& get_solution_names(void) const;
};


//static class
class Tecplot
{	
public:
	static void post_grid(const Grid<2>& grid) {
		This_::post_variables_.record_grid(grid);
		This_::file_writer_->write_grid_file(This_::post_variables_, This_::post_file_path_);
	}
	static void post_solution(const Discretized_Solution& discretized_solution) {
		const auto& solution_names = discretized_solution.get_solution_names();
		const auto post_point_solutions_by_variable = discretized_solution.calculate_post_point_solutions_by_variables();
		const auto num_solution_variable = solution_names.size();
		
		for (ushort i = 0; i < num_solution_variable; ++i) 
			This_::post_variables_.record_variable(solution_names[i], post_point_solutions_by_variable[i]);		
	}
	static void record_variables(const std::string_view name, const std::vector<double>& values) {
		This_::post_variables_.record_variable(name, values);
	}
		

private:
	Tecplot(void) = delete;

private:
	using This_ = Tecplot;

	static inline std::string post_file_path_;
	static inline Post_Variables post_variables_;
	static inline std::unique_ptr<Tecplot_File_Writer> file_writer_;
};































enum class Post_File_Type {
	grid,
	solution
};


enum class Post_File_Format {
	ASCII,
	binary
};






class Tecplot
{
private:
	Tecplot(void) = delete;

private:
	using This_ = Tecplot;
	
public:
	static inline bool post_condition_ = false;

private:
	static inline Post_File_Format file_format_ = Post_File_Format::binary;
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

public:
	static void set_path(const std::string& path) { This_::path_ = path; };
	static void syncronize_time(const double& current_time) { This_::time_ptr_ = &current_time; };

public:
	template <typename T>
	static void record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);

	static void record_cell_indexes(void);

	template <typename T>
	static void conditionally_record_cell_variables(const std::string& data_name, const std::vector<T>& cell_datas);

	static void conditionally_record_cell_indexes(void);

public:
	template <typename Governing_Equation>
	static void initialize(const ushort post_order, const Post_File_Format format);

	template <ushort space_dimension>
	static void post_grid(const Grid<space_dimension>& grid);

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
	static void write_ASCII_header(const Post_File_Type file_type, const std::string_view post_file_path);
	static void write_ASCII_grid_post_file(const std::vector<std::vector<double>>& post_coordinate_blocks, const std::vector<std::vector<int>>& connectivities);
	static void write_ASCII_solution_post_file(const std::vector<std::vector<double>>& post_solution_datas, const std::string& comment = "");


	static void write_binary_header(const Post_File_Type file_type, const std::string_view post_file_path);
	static void write_binary_grid_post_file(const std::vector<std::vector<double>>& post_coordinate_blocks, const std::vector<std::vector<int>>& connectivities);
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

namespace ms {
	std::string to_string(const Zone_Type zone_type);
}

//template definition part
template <typename Governing_Equation>
void Tecplot::initialize(const ushort post_order, const Post_File_Format format) {
	This_::file_format_ = format;
	This_::post_order_ = post_order;

	if constexpr (Governing_Equation::space_dimension() == 2) {
		if constexpr (ms::is_SCL<Governing_Equation>) {
			This_::grid_variables_str_ = "X,Y";
			This_::solution_variables_str_ = "q";
			This_::zone_type_ = Zone_Type::FETriangle;
		}
		else if constexpr (ms::is_Euler<Governing_Equation>) {
			This_::grid_variables_str_ = "X,Y";
			This_::solution_variables_str_ = "rho,rhou,rhov,rhoE,u,v,p";
			This_::zone_type_ = Zone_Type::FETriangle;
		}		
	}
	else if constexpr (Governing_Equation::space_dimension() == 3) {
		if constexpr (ms::is_SCL<Governing_Equation>) {
			This_::grid_variables_str_ = "X,Y,Z";
			This_::solution_variables_str_ = "q";
			This_::zone_type_ = Zone_Type::FETetrahedron;
		}
		else if constexpr (ms::is_Euler<Governing_Equation>) {
			This_::grid_variables_str_ = "X,Y,Z";
			This_::solution_variables_str_ = "rho,rhou,rhov,rhow,rhoE,u,v,w,p";
			This_::zone_type_ = Zone_Type::FETetrahedron;
		}		
	}
};

template <ushort space_dimension, typename Reconstruction_Method>
void Tecplot::initialize_HOM(const Grid<space_dimension>& grid, const Reconstruction_Method& reconstruct_method) {
	const auto set_of_post_nodes = grid.cell_set_of_post_nodes(This_::post_order_);
	
	const auto num_cell = set_of_post_nodes.size();
	This_::set_of_basis_post_points_.reserve(num_cell);

	for (uint i = 0; i < num_cell; ++i) 
		This_::set_of_basis_post_points_.push_back(reconstruct_method.basis_nodes(i, set_of_post_nodes[i]));
}

template <ushort space_dimension>
void Tecplot::post_grid(const Grid<space_dimension>& grid) {
	const auto set_of_post_nodes = grid.cell_set_of_post_nodes(This_::post_order_);
	
	const auto num_cell = set_of_post_nodes.size();
	This_::num_post_points_.reserve(num_cell);
	for (const auto& post_nodes : set_of_post_nodes)
		This_::num_post_points_.push_back(post_nodes.size());

	const auto post_coordinate_blocks = grid.cell_post_coordinate_blocks(set_of_post_nodes);
	const auto connectivities = grid.cell_set_of_connectivities(post_order_, set_of_post_nodes);

	This_::num_node_ = post_coordinate_blocks.front().size();
	This_::num_element_ = connectivities.size();

	if (This_::file_format_ == Post_File_Format::binary)
		This_::write_binary_grid_post_file(post_coordinate_blocks, connectivities);
	else
		This_::write_ASCII_grid_post_file(post_coordinate_blocks, connectivities);
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

	if (!This_::additioinal_data_name_to_values_.contains(data_name))
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

	if (This_::file_format_ == Post_File_Format::binary)
		This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
	else
		This_::write_ASCII_solution_post_file(post_solution_binary_datas, comment);	
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

	if (This_::file_format_ == Post_File_Format::binary)
		This_::write_binary_solution_post_file(post_solution_binary_datas, comment);
	else
		This_::write_ASCII_solution_post_file(post_solution_binary_datas, comment);
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
		//후에 추가적인 processig 부분이 다양하게 필요한 경우
		//switch 역활을 하는 enum class 만들고 
		//initialize 부분에서 Governing Equation에 따라서 switch 변수를 할당한뒤
		//performance loss가 있지만 if 문으로 처리해야 Tecplot:: << 이 형식을 유지할 수 있다.

		constexpr ushort space_dimension = num_equation - 2; //temporary code

		convert_to_binary_data.resize(2 * num_equation - 1);

		for (auto& post_point_values : convert_to_binary_data)
			post_point_values.reserve(This_::num_node_);

		for (uint i = 0; i < This_::num_node_; ++i) {
			const auto& cvariable = post_point_solutions[i];
			const auto pvariable = Euler<space_dimension>::conservative_to_primitive(cvariable);

			for (size_t j = 0; j < num_equation; ++j)
				convert_to_binary_data[j].push_back(cvariable[j]);

			for (size_t j = 0; j < num_equation - 1; ++j)
				convert_to_binary_data[num_equation + j].push_back(pvariable[j]);
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
