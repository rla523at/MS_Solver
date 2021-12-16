#pragma once
#include "Text.h"

#include <map>

using ushort = unsigned short;

class Configuration
{
public:
	Configuration(const std::string_view file_path);

public://Query
	ushort space_dimension(void) const;

	const std::string& get_governing_equation(void) const;
	const std::string& get_initial_condition(void) const;

	const std::vector<std::string>& get_grid_file_paths(void) const;
	const std::string& get_grid_file_type(void) const;

	const std::string& get_spatial_discrete_scheme(void) const;
	ushort solution_degree(void) const;
	const std::string& get_reconstruction_scheme(void) const;

	const std::string& get_numerical_flux(void) const;

	const std::string& get_time_discrete_scheme(void) const;

	const std::string& get_time_step_method(void) const;
	double CFL_number(void) const;
	double constant_dt(void) const;

	const std::string& get_solve_end_controller_type(void) const;
	double solve_end_time(void) const;
	size_t solve_end_iter(void) const;

	const std::string& get_solve_post_controller_type(void) const;
	double solve_post_time_step(void) const;
	size_t solve_post_iter_unit(void) const;

	const std::string& get_post_base_path(void) const;
	const std::string& get_post_file_format(void) const;
	const std::string& get_post_point_location(void) const;
	ushort post_order(void) const;

	const std::string& get_error_type(void) const;

	const double* advection_speeds_ptr(void) const;
	const double* wave_lengths_ptr(void) const;
	const double* periodic_lengths_ptr(void) const;

	std::string configuration_str(const std::string& grid_file_name) const;
	std::string governing_equation_str(void) const;
	std::string grid_file_name(const std::string& grid_file_path) const;
	std::string initial_condition_str(void) const;
	std::string post_folder_path_str(const std::string& grid_file_name) const;

	

private:
	template <typename ValueType>	ValueType get(const std::map<Sentence, Sentence>& name_to_value, const std::string_view config_name) const
	{
		REQUIRE(name_to_value.contains(ms::get_upper_case(config_name.data())), std::string(config_name.data()) + " should be exist in configuration file");

		const auto& value_sentence = name_to_value.at(ms::get_upper_case(config_name.data()));
		return value_sentence.to_value<ValueType>();
	};
	ushort find_solution_degree(void) const;
	Text read_file(const std::string_view file_path) const;
	void set_value(const Text& config_text);
	std::string solve_end_condition_str(void) const;
	std::string time_step_str(void) const;

private:
	static constexpr ushort max_space_dimension = 3;
	std::string date_time_str_;

	//datas in configration.dat
	ushort space_dimension_ = 0;

	std::string governing_equation_;
	std::string initial_condition_;

	//std::string grid_file_path_;
	std::vector<std::string> grid_file_paths_;
	std::string grid_file_type_;

	std::string spatial_discrete_scheme_;	
	ushort solution_degree_ = 0;
	std::string reconstruction_scheme_;

	std::string numerical_flux_;
	std::string time_discrete_scheme_;

	std::string time_step_method_;
	double CFL_number_ = 0.0;
	double constant_dt_ = 0.0;

	std::string solve_end_controller_type_;
	double solve_end_time_ = 0.0;
	size_t solve_end_iter_ = 0;

	std::string solve_post_controller_type_;
	double solve_post_time_step_ = 0.0;
	size_t solve_post_iter_unit_ = 0;

	std::string post_base_path_;
	std::string post_file_format_;
	std::string post_point_location_;
	ushort post_order_ = 0;

	std::string error_type_;

	//linear advection option
	std::array<double, max_space_dimension> advection_speeds_ = { 0 };
	//sine wave option
	std::array<double, max_space_dimension> wave_lengths_ = { 0 };
	//exact solution option
	std::array<double, max_space_dimension> periodic_lengths_ = { 0 };
};

namespace ms
{
	std::string date_time_string(void);
}