#pragma once

#include <string>
#include <stdexcept>

#define LOCATION location_str(__FILE__,__FUNCTION__,__LINE__)
#define REQUIRE(requirement, message) require(requirement, message, LOCATION)
#define EXCEPT(message) require(false, message, LOCATION)

inline std::string location_str(const std::string& file_name, const std::string& function_name, const int num_line) {
	std::string location_str = file_name;

	location_str.erase(location_str.begin(), location_str.begin() + location_str.rfind("\\") + 1);

	location_str = "File\t\t: " + location_str + "\n";
	location_str += "Function\t: " + function_name + "\n";
	location_str += "Line\t\t: " + std::to_string(num_line) + "\n";

	return location_str;
}
inline void require(const bool requirement, const std::string_view message, const std::string& location_str) {
	if (!requirement) {
		std::string exception_str = "\n\n============================EXCEPTION============================\n";
		exception_str += location_str;
		exception_str = exception_str + "Message\t\t: " + message.data();
		exception_str += "\n============================EXCEPTION============================\n\n";
		throw std::runtime_error(exception_str);
	}
}
