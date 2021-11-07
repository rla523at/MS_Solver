#pragma once
#include "Discrete_Solution.h"

class Time_Step_Calculator
{
public:
	virtual double calculate(const Discrete_Solution& discrete_solution, const Governing_Equation& governing_equation) const abstract;
};

class CFL : public Time_Step_Calculator
{
public:
    CFL(const double cfl, const Grid& grid);

public://Query
    double calculate(const Discrete_Solution& discrete_solution, const Governing_Equation& governing_equation) const override;

private:
    double cfl_;
    std::vector<double> cell_volumes_;
    std::vector<std::vector<double>> cell_projected_volumes_;
};

class Constant_Dt : public Time_Step_Calculator
{
public:
    Constant_Dt(const double constant_dt) : constant_dt_(constant_dt) {};

public:
    double calculate(const Discrete_Solution& discrete_solution, const Governing_Equation& governing_equation) const override;

private:
    double constant_dt_;
};

class Time_Step_Calculator_Factory
{
public:
    static std::unique_ptr<Time_Step_Calculator> make(const Configuration& configuration, const Grid& grid)
    {
        const auto type_name = configuration.get("time_step_calculator_type");
        if (ms::contains_icase(type_name, "CFL"))
        {
            const auto cfl_number = configuration.get<double>("CFL_number");
            return std::make_unique<CFL>(cfl_number, grid);
        }
        else if (ms::contains_icase(type_name, "Constant", "dt"))
        {
            const auto constant_dt = configuration.get<double>("constant_dt");
            return std::make_unique<Constant_Dt>(constant_dt);
        }
    }
};

//#include <type_traits>
//
//template <double time_step_constant>
//class TSM { //time step method
//private:
//	TSM(void) = delete;
//
//public:
//	static constexpr double constant() { return time_step_constant; };
//}; 
//
//
//template <double time_step_constant>
//class CFL : public TSM<time_step_constant>
//{
//private:
//	CFL(void) = delete;
//
//public:
//	static std::string name(void) { return "CFL " + std::to_string(time_step_constant); };
//};
//
//
//template <double time_step_constant>
//class Constant_Dt : public TSM<time_step_constant>
//{
//private:
//	Constant_Dt(void) = delete;
//
//public:
//	static std::string name(void) { return "Constant " + std::to_string(time_step_constant) + " Time"; };
//};
//
//
//namespace ms {
//	template <typename T>
//	inline constexpr bool is_time_step_method = std::is_base_of_v<TSM<T::constant()>, T>;
//}
