#pragma once
#include "Text.h"

#include <map>
#include <string>

class Configuration
{
public:
	//Configuration(const std::string& file_name) {
	//	const auto config_text = this->read_File(file_name);
	//	this->set_Value(config_text);
	//};

	template <typename ValueType>
	ValueType get(const std::string_view config_name) const {
		REQUIRE(this->name_to_value_str_.contains(ms::upper_case(config_name.data())), "name should be exist in configuration file");		
		const auto& value_str = name_to_value_str_.at(ms::upper_case(config_name));
		return ms::string_to_value<ValueType>(value_str);
	};
	std::string get(const std::string_view config_name) const {
		REQUIRE(this->name_to_value_str_.contains(ms::upper_case(config_name.data())), "name should be exist in configuration file");
		return name_to_value_str_.at(ms::upper_case(config_name.data()));
	};

	//std::string Case_Description(void) const;

private:
	Text read_file(const std::string_view file_path) {
		Text config_text;
					
		config_text.read(file_path);

		const std::vector<char> trim_character_set = { ' ', ',', '=', '\t', '\n', '\r' };
		const std::string comment_tag = "//";

		for (auto& sentence : config_text) {
			const auto comment_pos = sentence.find_position(comment_tag);
			sentence.remove_all_from_here(comment_pos);
		}

		config_text.remove_empty_line();
		return config_text;
	};

	void set_Value(const Text& config_text) {
		const auto delimiter = '=';
		
		for (const auto& sentence : config_text) {
			auto parsed_sentences = sentence.parse(delimiter);

			auto name_str = parsed_sentences.front().upper_case().to_string();
			auto value_sentence = parsed_sentences.back().to_string();

			this->name_to_value_str_.emplace(std::move(name_str), std::move(value_sentence));
		}
	}

private:
	std::map<std::string, std::string> name_to_value_str_;
};
