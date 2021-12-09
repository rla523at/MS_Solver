#pragma once

#include "Exception.h"

#include <algorithm>
#include <fstream>
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sstream>
#include <vector>

class Sentence
{
public:
	Sentence(void) = default;
	Sentence(const char* cp)
		:contents_(cp) {};
	Sentence(const std::string& str)
		:contents_(str) {};
	
public://command	
	Sentence& operator<<(const std::string& str);
	
	template<typename T>	Sentence& insert_with_space(const std::vector<T>& values)
	{
		for (const auto value : values)
		{
			this->insert_with_space(value);
		}

		return *this;
	}
	template<typename T>	Sentence& insert_with_space(const T value) 
	{
		this->contents_ += " " + std::to_string(value);
		return *this;
	}
	template<>	Sentence& insert_with_space(const double value);
	void remove_after(const std::string_view target);
	void remove_from_here(const size_t position);
	void remove(const std::string_view target);
	void remove(const std::vector<char> targets);
	void upper_case(void);

public://Query
	bool operator<(const Sentence& other) const;
	bool operator==(const Sentence& other) const;

	bool contain_icase(const char* target) const;
	template <typename... Args>		bool contain_icase(const Args... args) const;
	size_t find_position(const std::string_view target) const;
	Sentence get_remove(const std::string_view target) const;
	std::vector<Sentence> parse(const char delimiter) const;
	std::string to_string(void) const;
	template <typename T>	T to_value(void) const;
	Sentence get_upper_case(void) const;

private:
	std::string contents_;
};

class Text
{
public:
	Text(void) = default;
	Text(std::initializer_list<Sentence> list);

public://command
	Sentence& operator[](const size_t index);
	Text& operator<<(const std::string& str);
	Text& operator<<(std::string&& str);	

	void add_empty_lines(const size_t num_line);
	std::vector<Sentence>::iterator begin(void);
	void clear(void);
	std::vector<Sentence>::iterator end(void);
	void merge(Text&& other);
	void remove_empty_line(void);
	void read(const std::string_view read_file_path);
	void read(std::ifstream& file, const size_t num_read_line);

public://Query
	bool operator==(const Text& other) const;	
	void add_write(const std::string_view write_file_path) const;
	std::vector<Sentence>::const_iterator begin(void) const;
	std::vector<Sentence>::const_iterator end(void) const;
	void write(const std::string_view write_file_path) const;
	std::string to_string(void) const;
	size_t size(void) const;

private:
	std::vector<Sentence> senteces_;
};

class Binary_Writer
{
public:
	Binary_Writer(const std::string_view file_path);
	Binary_Writer(const std::string_view file_path, std::ios_base::openmode mode);

public://Command
	template <typename T>	Binary_Writer& operator<<(const T value);
	template <typename T>	Binary_Writer& operator<<(const std::vector<T>& values);
	template <>		Binary_Writer& operator<<(const char* value);
	template <>		Binary_Writer& operator<<(const std::string& str);

private:
	std::ofstream binary_file_stream_;
};

std::ostream& operator<<(std::ostream& ostream, const Sentence& sentece);
std::ostream& operator<<(std::ostream& ostream, const Text& text);
template <typename T>	Binary_Writer& Binary_Writer::operator<<(const T value) 
{
	this->binary_file_stream_.write(reinterpret_cast<const char*>(&value), sizeof(T));
	return *this;
}
template <typename T>	Binary_Writer& Binary_Writer::operator<<(const std::vector<T>& values) 
{
	for (const auto value : values)
		this->binary_file_stream_.write(reinterpret_cast<const char*>(&value), sizeof(T));
	return *this;
}

namespace ms
{
	bool contains_icase(const std::string& str, const char* target);
	bool contains_icase(const std::string& str, const std::string& target);
	template <typename... Args>		bool contains_icase(const std::string& str, const Args... args)
	{
		//static_assert((... && std::is_same_v<Args, const char*>), "every arguments should be array of char");
		return (ms::contains_icase(str, args) && ...);
	};
	std::string double_to_string(const double val);
	std::string double_to_str_sp(const double value); //double to string with show point
	std::vector<std::string> file_paths_in_path(const std::string& path);
	size_t find_icase(const std::string& str, const std::string& target);
	std::string get_replace(const std::string& str, const std::string_view target, const std::string_view replacement);
	std::string get_remove(const std::string& str, const char target);
	std::string get_remove(const std::string& str, const std::string_view target);
	std::string get_upper_case(const std::string& str);
	void make_path(std::string_view file_path);
	std::vector<std::string> parse(const std::string& str, const char delimiter);
	std::vector<std::string> parse(const std::string& str, const std::vector<char>& delimiters);
	void replace(std::string& str, const char target, const char replacement);
	void replace(std::string& str, const std::string_view target, const std::string_view replacement);
	void remove(std::string& str, const char target);
	void remove(std::string& str, const std::string_view target);
	void rename(const std::string& path, const std::string& old_name, const std::string& new_name);
	size_t rfind_nth(const std::string& object_str, const std::string& target_str, const size_t n);
	template<typename T>	T string_to_value(const std::string& str)
	{
		std::istringstream iss(str);
		T value;
		iss >> value;
		return value;
	};
	template<typename T>	std::vector<T> string_to_value_set(const std::vector<std::string>& str_set)
	{
		std::vector<T> result;
		result.reserve(str_set.size());

		for (const auto& str : str_set)
			result.push_back(ms::string_to_value<T>(str));

		return result;
	};
	template<typename T>	std::vector<T> sentences_to_values(const std::vector<Sentence>& sentences)
	{
		std::vector<T> result;
		result.reserve(sentences.size());

		for (const auto& sentence : sentences)
			result.push_back(sentence.to_value<T>());

		return result;
	};
	void upper_case(std::string& str);
}


//template definition part
template <typename... Args>		bool Sentence::contain_icase(const Args... args) const
{
	return ms::contains_icase(this->contents_, args...);
};
template <typename T>	T Sentence::to_value(void) const
{
	return ms::string_to_value<T>(contents_);
}

