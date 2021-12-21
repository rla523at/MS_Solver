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