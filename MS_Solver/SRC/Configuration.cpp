#include "../INC/Configuration.h"

Configuration::Configuration(const std::string_view file_path) 
{
	const auto config_text = this->read_file(file_path);
	this->set_value(config_text);
};

std::string Configuration::get(const std::string_view config_name) const
{
	REQUIRE(this->name_to_value_str_.contains(ms::get_upper_case(config_name.data())), "name should be exist in configuration file");
	return name_to_value_str_.at(ms::get_upper_case(config_name.data()));
};

//std::string Case_Description(void) const;


Text Configuration::read_file(const std::string_view file_path) const
{
	Text config_text;
	config_text.read(file_path);

	const std::vector<char> trim_character_set = { ' ', '\t' };
	const std::string comment_tag = "//";

	for (auto& sentence : config_text) 
	{
		sentence.remove(trim_character_set);

		const auto comment_pos = sentence.find_position(comment_tag);

		if (comment_pos != std::string::npos)
		{
			sentence.remove_from_here(comment_pos);
		}
	}

	config_text.remove_empty_line();
	return config_text;
};

void Configuration::set_value(const Text& config_text) 
{
	const auto delimiter = '=';

	for (const auto& sentence : config_text) 
	{
		auto parsed_sentences = sentence.parse(delimiter);

		auto name_str = parsed_sentences.front().get_upper_case().get_string();
		auto value_str = parsed_sentences.back().get_string();

		this->name_to_value_str_.emplace(std::move(name_str), std::move(value_str));
	}
}