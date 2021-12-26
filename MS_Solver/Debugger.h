#pragma once
#include <array>
#include <iostream>
#include <vector>

class Debugger
{
public:
	inline static size_t count_ = 0;
	inline static std::array<bool, 10> conditions_ = { false };

private:
	Debugger(void) = delete;

private:
	using This_ = Debugger;
};


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec)
{
	for (const auto& elem : vec)
		os << elem << " ";
	return os;
}