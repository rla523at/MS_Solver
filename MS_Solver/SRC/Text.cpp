#include "../INC/Text.h"

Text::Text(std::ifstream& file, const size_t num_read_line) {
	this->read(file, num_read_line);
}

Text& Text::operator<<(const std::string& str) {
	this->push_back(str);
	return *this;
}

Text& Text::operator<<(std::string&& str) {
	this->push_back(std::move(str));
	return *this;
}

void Text::add_write(const std::string& file_path) const {
	ms::make_path(file_path);
	std::ofstream output_file(file_path, std::ios::app);
	dynamic_require(output_file.is_open(), "Fail to open file" + file_path);

	const auto num_sentence = this->size();
	for (auto i = this->begin(); i != this->end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->back();

	output_file.close();
}

Text& Text::read_line_by_line(const std::string& file_path) {	
	std::ifstream file_stream(file_path);
	dynamic_require(file_stream.is_open(), "Fail to open file" + file_path);

	std::string str;
	while (std::getline(file_stream, str))
		this->push_back(std::move(str));

	file_stream.close();
	return *this;
}

void Text::read(std::ifstream& file_stream, const size_t num_read_line) {
	dynamic_require(file_stream.is_open(), "Fail to open file");

	size_t index = 0;
	std::string str;
	while (std::getline(file_stream, str)) {
		this->push_back(std::move(str));
		if (++index == num_read_line)
			break;
	}

	//file_stream.close();
}

Text& Text::remove_empty_line(void) {
	this->erase(std::remove(this->begin(), this->end(), ""), this->end());
	return *this;
}

void Text::write(const std::string& file_path) const {
	ms::make_path(file_path);
	std::ofstream output_file(file_path);
	dynamic_require(output_file.is_open(), "Fail to open file" + file_path); 

	const auto num_sentence = this->size();
	for (auto i = this->begin(); i != this->end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->back();

	output_file.close();
}

void Text::merge(Text&& other) {
	this->insert(this->end(), std::make_move_iterator(other.begin()), std::make_move_iterator(other.end()));
}

std::ostream& operator<<(std::ostream& ostream, const Text& text) {
	for (const auto& line : text)
		ostream << line << "\n";
	return ostream;
}

Binary_Writer::Binary_Writer(const std::string_view file_path) {
	ms::make_path(file_path);
	binary_file_stream_.open(file_path.data(), std::ios::binary | std::ios::app);

	dynamic_require(this->binary_file_stream_.is_open(), "file should be opened");
}

template <>
Binary_Writer& Binary_Writer::operator<<(const char* value) {
	this->binary_file_stream_.write(reinterpret_cast<const char*>(value), sizeof(value));
	return *this;
}

template <>
Binary_Writer& Binary_Writer::operator<<(const std::string& str) {
	this->binary_file_stream_ << str;
	return *this;
}

namespace ms {
	std::vector<std::string> parse(const std::string& str, const char delimiter) {
		if (str.empty())
			return std::vector<std::string>();

		std::vector<std::string> parsed_string_set;
		for (auto iter1 = str.begin();;) {
			auto iter2 = std::find(iter1, str.end(), delimiter);

			if (iter1 != iter2)
				parsed_string_set.emplace_back(iter1, iter2);

			if (iter2 == str.end())
				return parsed_string_set;

			iter1 = iter2 + 1;
		}
	};

	std::string erase(const std::string& str, const std::string& target) {
		const auto target_size = target.size();

		auto result = str;
		while (true) {
			const auto position = result.find(target);

			if (position == std::string::npos)
				return result;

			result.erase(position, target_size);
		}
	}

	size_t find_icase(const std::string& str, const std::string& target) {
		auto u_str = ms::upper_case(str);
		auto u_target = ms::upper_case(target);

		return u_str.find(u_target);
	}

	std::string upper_case(const std::string& str) {
		auto result = str;
		std::transform(result.begin(), result.end(), result.begin(), toupper);
		//std::transform(result.begin(), result.end(), result.begin(), std::toupper); //filesystem이랑 있으면 충돌
		return result;
	}

	bool is_there_icase(const std::string& str, const std::string& target) {
		return ms::find_icase(str, target) != std::string::npos;
	}

	std::string double_to_str_sp(const double value) {
		std::ostringstream os;
		os << std::setprecision(16) << std::showpoint << value;
		return os.str();
	}

	Text extract_file_path_text(const std::string& path) {
		Text file_name_text;

		std::filesystem::directory_iterator iter(path);
		while (iter != std::filesystem::end(iter)) {
			const auto& entry = *iter;

			if (entry.is_directory()) {
				iter++;
				continue;
			}

			file_name_text << entry.path().string();
			iter++;
		}

		return file_name_text;
	}

	void make_path(std::string_view file_path) {
		const auto file_name_size = file_path.size() - file_path.find_last_of("/") - 1;
		file_path.remove_suffix(file_name_size);

		if (file_path.empty())
			return;

		std::filesystem::path p(file_path);
		if (std::filesystem::exists(p))
			return;
		else
			std::filesystem::create_directories(p);
	}
}


