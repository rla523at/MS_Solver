#pragma once
#include "Text.h"

#include <iostream>
#include <sstream>


class Log
{
private:
	inline static std::string path_;
	inline static Text log_txt_;

public:
	inline static std::ostringstream content_;

private:
	Log(void) = delete;

public:
	static void set_path(const std::string& path) {	
		Log::path_ = path;
	}

	static void print(void) {
		auto log_str = Log::content_.str();

		std::cout << log_str;
		Log::log_txt_ << std::move(log_str);

		std::ostringstream tmp;
		Log::content_.swap(tmp);
	}

	static void write_error_text(const std::string& error_string) {		
		const auto pos = ms::rfind_nth(Log::path_, "/", 2);
		const auto sub_str = Log::path_.substr(pos + 1);
		const auto parsed_str = ms::parse(sub_str, '_');
		const auto grid_file_name = parsed_str[0];

		Text txt = { grid_file_name + "    \t" + error_string };

		auto error_text_path = Log::path_;
		error_text_path.erase(pos + 1);

		txt.add_write(error_text_path);
	}

	static void write(void) {
		const auto file_path = Log::path_ + "_log.txt";
		Log::log_txt_.write(file_path);
		Log::log_txt_.clear();	
	}

	static std::string date_string(void) {
		time_t date_calculator = time(nullptr);
		tm date;
		localtime_s(&date, &date_calculator);
		
		return std::to_string(date.tm_year + 1900) + "." + std::to_string(date.tm_mon + 1) + "." + std::to_string(date.tm_mday) 
			+ "__(" + std::to_string(date.tm_hour) + "£º" + std::to_string(date.tm_min) + ")";
	}
};

