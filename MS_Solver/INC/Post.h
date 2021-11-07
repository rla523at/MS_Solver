#pragma once
#include "Discretized_Solution.h"
#include "Tecplot.h"

//static class
class Post_Processing
{
public:
	static void initialize(const Configuration& configuration);
	static void post_grid(const Grid& grid);
	static void post_solution(const Discrete_Solution& discretized_solution);
	static void record_variables(const std::string_view name, const std::vector<double>& values);
	static void syncronize_solution_time(const double& get_solution_time);

private:
	Post_Processing(void) = delete;

private:
	using This_ = Post_Processing;

	static inline std::string post_file_path_;
	static inline ushort post_order_ = 0;
	static inline Post_Variables post_variables_;
	static inline Tecplot_File_Writer file_writer_;
};