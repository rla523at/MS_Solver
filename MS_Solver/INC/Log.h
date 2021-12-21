#pragma once
#include "Configuration.h"

#include <iostream>
#include <sstream>

#define LOG Log::get_instance()

class Log
{
private:
	enum class Command
	{
		print,
		print_off,
		print_on
	};

public:
	static Log& get_instance(void);

public://Command
	template <typename T> Log& operator<<(const T& value)
	{
		this->content_ << value;
		return *this;
	}
	Log& operator<<(const Command& command);
	
	void clear(void);	
	void set_configuration(const Configuration& configuration);
	void set_path(const std::string& grid_folder_path);

public://Query
	void write(void) const;
	void write_error_text(const std::string& grid_file_name, const std::vector<double>& error_values) const;

private:
	void print(void);
	std::string make_error_file_path(const std::string& path);

private:
	bool do_write = true;

	std::ostringstream content_;
	bool is_print_on = true;
	Text log_txt_;
	std::string log_folder_path_;
	std::string error_log_folder_path_;

public:
	static constexpr Command print_ = Command::print;
	static constexpr Command print_off = Command::print_off;
	static constexpr Command print_on = Command::print_on;
};

