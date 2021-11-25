#include "../INC/Text.h"

Sentence& Sentence::operator<<(const std::string& str) 
{
	this->contents_ += str;
	return *this;
}

template<>
Sentence& Sentence::insert_with_space(const double value) 
{
	this->contents_ += " " + ms::double_to_string(value);
	return *this;
}

void Sentence::remove_after(const std::string_view target) 
{
	const auto pos = this->contents_.find(target.data());
	this->contents_.erase(pos + 1);
}

void Sentence::remove_from_here(const size_t position) 
{
	REQUIRE(position < this->contents_.size(), "position can not exceed given range");
	this->contents_.erase(position);
}

void Sentence::remove(const std::string_view target)
{
	ms::remove(this->contents_, target);
}

void Sentence::remove(const std::vector<char> targets) 
{
	for (const auto target : targets)
		ms::remove(this->contents_, target);
}

void Sentence::upper_case(void) 
{
	ms::upper_case(this->contents_);
}

bool Sentence::operator==(const Sentence& other) const 
{
	return this->contents_ == other.contents_;
}

bool Sentence::contain_icase(const char* target) const
{
	return ms::contains_icase(this->contents_, target);
}

size_t Sentence::find_position(const std::string_view target) const 
{
	return this->contents_.find(target.data());
}

Sentence Sentence::get_remove(const std::string_view target) const
{
	auto result = *this;
	result.remove(target);
	return result;
}


std::vector<Sentence> Sentence::parse(const char delimiter) const 
{
	auto parsed_strs =ms::parse(this->contents_, delimiter);
	
	std::vector<Sentence> parsed_senteces;
	parsed_senteces.reserve(parsed_strs.size());

	for (auto& str : parsed_strs)
		parsed_senteces.push_back(std::move(str));

	return parsed_senteces;
}

std::string Sentence::get_string(void) const 
{
	return this->contents_;
}

Sentence Sentence::get_upper_case(void) const 
{	
	return ms::get_upper_case(this->contents_);
}

Text::Text(std::initializer_list<std::string> list) 
{
	this->senteces_.reserve(list.size());

	for (const auto& str : list)
		this->senteces_.push_back(str);
}

Sentence& Text::operator[](const size_t index) 
{
	REQUIRE(index < this->size(), "index can not exceed given range");
	return this->senteces_[index];
}

Text& Text::operator<<(const std::string& str) 
{
	this->senteces_.push_back(str);
	return *this;
}

Text& Text::operator<<(std::string&& str) 
{
	this->senteces_.push_back(std::move(str));
	return *this;
}

void Text::add_empty_lines(const size_t num_line) 
{
	this->senteces_.resize(this->senteces_.size() + num_line);
}

std::vector<Sentence>::iterator Text::begin(void) 
{
	return this->senteces_.begin();
}

void Text::clear(void)
{
	this->senteces_.clear();
}


std::vector<Sentence>::iterator Text::end(void) 
{
	return this->senteces_.end();
}

void Text::merge(Text&& other) 
{
	this->senteces_.insert(this->senteces_.end(), std::make_move_iterator(other.senteces_.begin()), std::make_move_iterator(other.senteces_.end()));
}

void Text::read(const std::string_view file_path) 
{
	std::ifstream file_stream(file_path);
	REQUIRE(file_stream.is_open(), "Fail to open file");

	std::string str;
	while (std::getline(file_stream, str))
		this->senteces_.push_back(std::move(str));

	file_stream.close();
}

void Text::read(std::ifstream& file_stream, const size_t num_read_line) 
{
	REQUIRE(file_stream.is_open(), "Fail to open file");

	size_t index = 0;
	std::string str;
	while (std::getline(file_stream, str)) 
{
		this->senteces_.push_back(std::move(str));
		if (++index == num_read_line)
			break;
	}
}

void Text::remove_empty_line(void) 
{
	this->senteces_.erase(std::remove(this->senteces_.begin(), this->senteces_.end(), ""), this->senteces_.end());
}
bool Text::operator==(const Text& other) const 
{
	return this->senteces_ == other.senteces_;
}
void Text::add_write(const std::string_view file_path) const 
{
	ms::make_path(file_path);
	std::ofstream output_file(file_path, std::ios::app);
	REQUIRE(output_file.is_open(), "output file stream should be opend before write");

	for (auto i = this->senteces_.begin(); i != this->senteces_.end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->senteces_.back();

	output_file.close();
}
std::vector<Sentence>::const_iterator Text::begin(void) const 
{
	return this->senteces_.begin();
}
std::vector<Sentence>::const_iterator Text::end(void) const 
{
	return this->senteces_.end();
}
void Text::write(const std::string_view file_path) const 
{
	ms::make_path(file_path);
	std::ofstream output_file(file_path);
	REQUIRE(output_file.is_open(), "output file stream should be opend before write");

	for (auto i = this->senteces_.begin(); i != this->senteces_.end() - 1; ++i)
		output_file << *i << "\n";
	output_file << this->senteces_.back();

	output_file.close();
}
std::string Text::to_string(void) const 
{
	std::string str;
	for (auto i = this->senteces_.begin(); i != this->senteces_.end() - 1; ++i)
		str += i->get_string() + "\n";
	str += this->senteces_.back().get_string();

	return str;
}
size_t Text::size(void) const 
{
	return this->senteces_.size();
}


Binary_Writer::Binary_Writer(const std::string_view file_path) 
{
	ms::make_path(file_path);
	binary_file_stream_.open(file_path.data(), std::ios::binary);

	REQUIRE(this->binary_file_stream_.is_open(), "file should be opened");
}

Binary_Writer::Binary_Writer(const std::string_view file_path, std::ios_base::openmode mode) 
{
	ms::make_path(file_path);
	binary_file_stream_.open(file_path.data(), std::ios::binary | mode);

	REQUIRE(this->binary_file_stream_.is_open(), "file should be opened");
}

template <>
Binary_Writer& Binary_Writer::operator<<(const char* value) 
{
	this->binary_file_stream_.write(reinterpret_cast<const char*>(value), sizeof(value));
	return *this;
}

template <>
Binary_Writer& Binary_Writer::operator<<(const std::string& str) 
{
	this->binary_file_stream_ << str;
	return *this;
}

namespace ms 
{
	bool contains_icase(const std::string& str, const char* target) 
	{
		return ms::find_icase(str, target) != std::string::npos;
	}

	bool contains_icase(const std::string& str, const std::string& target)
	{
		return ms::find_icase(str, target) != std::string::npos;
	}

	std::string double_to_string(const double val) 
	{
		constexpr size_t precision = 16;
		std::stringstream stream;
		stream << std::setprecision(precision) << std::noshowpoint << val;
		return stream.str();
	}

	std::string double_to_str_sp(const double value) 
	{
		std::ostringstream os;
		os << std::setprecision(16) << std::showpoint << value;
		return os.str();
	}

	std::vector<std::string> file_paths_in_path(const std::string& path) 
	{
		std::vector<std::string> file_name_text;

		std::filesystem::directory_iterator iter(path);
		while (iter != std::filesystem::end(iter)) 
		{
			const auto& entry = *iter;

			if (entry.is_directory()) 
			{
				iter++;
				continue;
			}

			file_name_text.push_back(entry.path().string());
			iter++;
		}

		return file_name_text;
	}

	size_t find_icase(const std::string& str, const std::string& target) 
	{
		auto u_str = ms::get_upper_case(str);
		auto u_target = ms::get_upper_case(target);

		return u_str.find(u_target);
	}

	std::string get_replace(const std::string& str, const std::string_view target, const std::string_view replacement)
	{
		auto result = str;
		ms::replace(result, target, replacement);
		return result;
	}

	std::string get_remove(const std::string& str, const char target)
	{
		auto result = str;
		ms::remove(result, target);
		return result;
	}


	std::string get_remove(const std::string& str, const std::string_view target)
	{
		auto result = str;
		ms::remove(result, target);
		return result;
	}

	std::string get_upper_case(const std::string& str) 
	{
		auto result = str;
		ms::upper_case(result);
		return result;
	}

	void make_path(std::string_view file_path)
	{
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

	std::vector<std::string> parse(const std::string& str, const char delimiter) 
	{
		if (str.empty())
			return std::vector<std::string>();

		std::vector<std::string> parsed_string_set;
		for (auto iter1 = str.begin();;) 
		{
			auto iter2 = std::find(iter1, str.end(), delimiter);

			if (iter1 != iter2)
				parsed_string_set.emplace_back(iter1, iter2);

			if (iter2 == str.end())
				return parsed_string_set;

			iter1 = iter2 + 1;
		}
	};

	std::vector<std::string> parse(const std::string& str, const std::vector<char>& delimiters) 
	{
		if (delimiters.empty())
			return 
{ str };

		const auto num_delimiter = delimiters.size();

		const auto reference_delimiter = delimiters[0];

		auto temp_str = str;
		for (size_t i = 1; i < num_delimiter; ++i)
		{
			ms::replace(temp_str, delimiters[i], reference_delimiter);
		}

		return ms::parse(temp_str, reference_delimiter);
	}

	void replace(std::string& str, const char target, const char replacement) 
	{
		while (true) 
		{
			const auto pos = str.find(target);

			if (pos == std::string::npos)
				break;

			str[pos] = replacement;
		}
	}

	void replace(std::string& str, const std::string_view target, const std::string_view replacement)
	{
		if (target.empty())
			return;

		while (true) 
		{
			const auto pos = str.find(target.data());

			if (pos == std::string::npos)
				break;

			str.replace(pos, target.size(), replacement.data());
		}
	}

	void remove(std::string& str, const char target) 
	{
		while (true) 
		{
			const auto pos = str.find(target);

			if (pos == std::string::npos)
				break;

			str.erase(pos, 1);
		}
	}

	void remove(std::string& str, const std::string_view target) 
	{
		ms::replace(str, target, "");
	}

	void rename(const std::string& path, const std::string& old_name, const std::string& new_name)
	{
		std::filesystem::rename(path + old_name, path + new_name);
	}

	size_t rfind_nth(const std::string& object_str, const std::string& target_str, const size_t n) 
{
		if (n < 1)
			return std::string::npos;

		auto pos = std::string::npos;
		for (size_t i = 0; i < n; ++i) 
{
			if (pos == 0)
				return std::string::npos;

			pos = object_str.rfind(target_str, pos - 1); //rfind next pos
		}

		return pos;
	}

	void upper_case(std::string& str) 
	{
		std::transform(str.begin(), str.end(), str.begin(), toupper);
		//std::transform(result.begin(), result.end(), result.begin(), std::toupper); //filesystem이랑 있으면 충돌
	}
}



std::ostream& operator<<(std::ostream& ostream, const Sentence& sentece) 
{
	return ostream << sentece.get_string();
}

std::ostream& operator<<(std::ostream& ostream, const Text& text) 
{
	return ostream << text.to_string();
}