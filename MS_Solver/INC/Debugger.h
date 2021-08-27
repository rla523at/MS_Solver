#pragma once
#include <vector>
#include <iostream>

class Debugger
{
private:
	using This_ = Debugger;

public:
	inline static size_t count_ = 0;
	inline static std::vector<bool> conditions_;

private:
	Debugger(void) = delete;
};


template <typename T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& vec) {
	for (const auto& elem : vec)
		os << elem << " ";
	return os;
}