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
			+ "_" + std::to_string(date.tm_hour) + "h" + std::to_string(date.tm_min) + "m";		
	}
};

