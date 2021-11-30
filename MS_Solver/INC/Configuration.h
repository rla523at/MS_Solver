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
		REQUIRE(this->name_to_value_str_.contains(ms::get_upper_case(config_name.data())), "name should be exist in configuration file");		
		const auto& value_str = name_to_value_str_.at(ms::get_upper_case(config_name.data()));
		return ms::string_to_value<ValueType>(value_str);
	};

	//std::string Case_Description(void) const;

private:
	Text read_file(const std::string_view file_path) const;
	void set_value(const Text& config_text);

private:
	std::map<std::string, std::string> name_to_value_str_;
};
