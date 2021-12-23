#pragma once
#include "Post_Variable_Convertor.h"
#include "Tecplot_File_Writer.h"

class Post_Processor//static class
{
public:
	static void set_configuration(const Configuration& configuration);
	static void initialize(const Configuration& configuration, const Grid& grid, Discrete_Solution& discrete_solution);
	static void post_grid(void);
	static void post_solution(const std::string_view sv = "");
	static void record_solution(void);
	static void record_variables(const std::string& name, const std::vector<double>& values);	
	static void set_path(const std::string& path);
	static void syncronize_solution_time(const double& solution_time);

private:
	Post_Processor(void) = delete;

private:
	using This_ = Post_Processor;
	
	static inline bool do_write_ = true;
	static inline bool is_initialized_ = false;
	static inline bool is_solution_recorded_ = false;
	static inline std::string post_folder_path_;
	static inline ushort post_order_ = 0;
	static inline std::unique_ptr<Post_Variables> post_variables_;
	static inline std::unique_ptr<Tecplot_File_Writer> file_writer_;
};	
