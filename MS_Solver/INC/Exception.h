#pragma once
#include <string>
#include <stdexcept>
#include <sstream>

//#define DEVELOPE

#ifdef DEVELOPE
#define LOCATION ms::location_str(__FILE__,__FUNCTION__,__LINE__)
#define REQUIRE(requirement, message) ms::require(requirement, message, LOCATION)
#define EXCEPTION(message) ms::require(false, message, LOCATION)
#else
#define LOCATION 
#define REQUIRE(requirement, message)
#define EXCEPTION(message)
#endif // DEVELOPE





namespace ms
{
	inline std::string location_str(const std::string& file_name, const std::string& function_name, const int num_line)
	{
		std::string location_str = file_name;

		location_str.erase(location_str.begin(), location_str.begin() + location_str.rfind("\\") + 1);

		location_str = "File\t\t: " + location_str + "\n";
		location_str += "Function\t: " + function_name + "\n";
		location_str += "Line\t\t: " + std::to_string(num_line) + "\n";

		return location_str;
	}
	inline void require(const bool requirement, const std::string_view message, const std::string& location_str = "")
	{
		if (!requirement)
		{
			std::ostringstream os;

			os << "==============================EXCEPTION========================================\n";
			os << location_str;
			os << "Message\t\t: " << message.data() << "\n";
			os << "==============================EXCEPTION========================================\n\n";

			throw std::runtime_error(os.str());
		}
	}
}