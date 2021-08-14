#pragma once
#include <type_traits>

template <double time_step_constant>
class TSM { //time step method
private:
	TSM(void) = delete;

public:
	static constexpr double constant() { return time_step_constant; };
}; 


template <double time_step_constant>
class CFL : public TSM<time_step_constant>
{
private:
	CFL(void) = delete;

public:
	static std::string name(void) { return "CFL " + std::to_string(time_step_constant); };
};


template <double time_step_constant>
class Constant_Dt : public TSM<time_step_constant>
{
private:
	Constant_Dt(void) = delete;

public:
	static std::string name(void) { return "Constant " + std::to_string(time_step_constant) + " Time"; };
};


namespace ms {
	template <typename T>
	inline constexpr bool is_time_step_method = std::is_base_of_v<TSM<T::constant()>, T>;
}
