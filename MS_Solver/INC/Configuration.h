#pragma once
#include "Text.h"

#include <map>

class Configuration
{
public:
	Configuration(const std::string_view file_path);

public://Query
	std::string get(const std::string_view config_name) const;
	template <typename ValueType>	ValueType get(const std::string_view config_name) const
	{
		REQUIRE(this->name_to_value_str_.contains(ms::get_upper_case(config_name.data())), std::string(config_name.data()) + " should be exist in configuration file");		
		const auto& value_str = name_to_value_str_.at(ms::get_upper_case(config_name.data()));
		return ms::string_to_value<ValueType>(value_str);
	};
	std::string post_folder_path_str(void) const;
	std::string configuration_str(void) const;

private:
	std::string governing_equation_str(void) const;
	std::string grid_file_name(void) const;
	std::string initial_condition_str(void) const;
	Text read_file(const std::string_view file_path) const;
	void set_value(const Text& config_text);
	std::string solve_end_condition_str(void) const;
	std::string time_step_str(void) const;

private:
	std::map<std::string, std::string> name_to_value_str_;
};

namespace ms
{
	std::string date_string(void);
}