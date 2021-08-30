#pragma once
#include <algorithm>
#include <fstream>	//file stream
#include <filesystem>
#include <stdexcept>
#include <string>
#include <string_view>
#include <sstream>
#include <vector>

#define dynamic_require(requirement, state) if (!(requirement)) throw std::runtime_error(state)

class Text : public std::vector<std::string>
{
public:
	template <typename ... Vals>
	explicit Text(Vals&&... values) : std::vector<std::string>(std::forward<Vals>(values)...) {};
	Text(std::initializer_list<std::string> list) : std::vector<std::string>( list ) {};
	Text(std::ifstream& file, const size_t num_read_line);

public:
	Text& operator<<(const std::string& str);
	Text& operator<<(std::string&& str);

public:
	void merge(Text&& other);
	Text& remove_empty_line(void);
	Text& read_line_by_line(const std::string& read_file_path);
	void read(std::ifstream& file, const size_t num_read_line);

public:
	void add_write(const std::string& write_file_path) const;
	void write(const std::string& write_file_path) const;
};

std::ostream& operator<<(std::ostream& ostream, const Text& text);


class Binary_Writer
{
private:
	std::ofstream binary_file_stream_;

public:
	Binary_Writer(const std::string_view file_path);
	Binary_Writer(const std::string_view file_path, std::ios_base::openmode mode);

	template <typename T>
	Binary_Writer& operator<<(const T value);

	template <typename T>
	Binary_Writer& operator<<(const std::vector<T>& values);

	template <>
	Binary_Writer& operator<<(const char* value);

	template <>
	Binary_Writer& operator<<(const std::string& str);
};


namespace ms {
	std::vector<std::string> parse(const std::string& str, const char delimiter);
	template<typename T>
	T string_to_value(const std::string& str);
	template<typename T>
	std::vector<T> string_to_value_set(const std::vector<std::string>& str_set);
	std::string remove(const std::string& str, const std::string& target);
	std::string upper_case(const std::string& str);
	size_t find_icase(const std::string& str, const std::string& target);
	size_t rfind_nth(const std::string& object_str, const std::string& target_str, const size_t n);
	bool is_there_icase(const std::string& str, const std::string& target);
	std::string double_to_str_sp(const double value); //double to string with show point
	Text extract_file_path_text(const std::string& path);
	void make_path(std::string_view file_path);
}


//template definition
template <typename T>
Binary_Writer& Binary_Writer::operator<<(const T value) {
	this->binary_file_stream_.write(reinterpret_cast<const char*>(&value), sizeof(T));
	return *this;
}

template <typename T>
Binary_Writer& Binary_Writer::operator<<(const std::vector<T>& values) {
	for (const auto value : values)
		this->binary_file_stream_.write(reinterpret_cast<const char*>(&value), sizeof(T));
	return *this;
}

namespace ms {
	template<typename T>
	T string_to_value(const std::string& str) {
		std::istringstream iss(str);
		T value;
		iss >> value;
		return value;
	};
	template<typename T>
	std::vector<T> string_to_value_set(const std::vector<std::string>& str_set) {
		std::vector<T> result;
		result.reserve(str_set.size());

		for (const auto& str : str_set)
			result.push_back(ms::string_to_value<T>(str));

		return result;
	};
}