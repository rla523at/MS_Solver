#pragma once
#include "Text.h"

#include <iostream>
#include <sstream>

#define LOG Log::get_instance()

class Log
{
private:
	enum class Command
	{
		print
	};

public://Command
	template <typename T> Log& operator<<(const T& value)
	{
		this->content_ << value;
		return *this;
	}
	Log& operator<<(const Command& command)
	{
		switch (command)
		{
		case Command::print:	this->print(); break;
		default:				EXCEPTION("not supported command");
		}
		return *this;
	}
	
	void clear(void)
	{
		this->log_txt_.clear();
	}
	static Log& get_instance(void)
	{
		static Log instance;
		return instance;
	}
	void set_path(const std::string& path) 
	{	
		this->path_ = path;
	}

public://Query
	void write(void) const
	{
		const auto file_path = this->path_ + "_log.txt";
		this->log_txt_.write(file_path);
	}
	void write_error_text(const std::string& error_string) const
	{
		const auto pos = ms::rfind_nth(this->path_, "/", 2);
		const auto sub_str = this->path_.substr(pos + 1);
		const auto parsed_str = ms::parse(sub_str, '_');
		const auto grid_file_name = parsed_str[0];

		Text txt = { grid_file_name + "    \t" + error_string };

		auto error_text_path = this->path_;
		error_text_path.erase(pos + 1);

		txt.add_write(error_text_path + "error.txt");
	}
	std::string date_string(void) const
	{
		time_t date_calculator = time(nullptr);
		tm date;
		localtime_s(&date, &date_calculator);

		return std::to_string(date.tm_year + 1900) + "." + std::to_string(date.tm_mon + 1) + "." + std::to_string(date.tm_mday)
			+ "__(" + std::to_string(date.tm_hour) + "£º" + std::to_string(date.tm_min) + ")";
	}

private:
	void print(void) 
	{
		auto log_str = this->content_.str();
		std::cout << log_str;

		this->log_txt_ << std::move(log_str);

		std::ostringstream tmp;
		this->content_.swap(tmp);
	}

private:
	std::string path_;
	Text log_txt_;
	std::ostringstream content_;

public:
	static constexpr Command print_ = Command::print;
};

