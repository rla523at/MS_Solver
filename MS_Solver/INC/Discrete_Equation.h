#pragma once
#include "Solve_Controller.h"
#include "Time_Discrete_Scheme.h"

class Discrete_Equation
{
public:
	Discrete_Equation(std::unique_ptr<Semi_Discrete_Equation>&& semi_discrete_equation, std::unique_ptr<Time_Discrete_Scheme>&& time_discrete_scheme,
		std::unique_ptr<Solve_End_Controller>&& end_controller, std::unique_ptr<Solve_Post_Controller>&& post_controller)
		: semi_discrete_equation_(std::move(semi_discrete_equation))
		, time_discrete_scheme_(std::move(time_discrete_scheme))
		, end_controller_(std::move(end_controller))
		, post_controller_(std::move(post_controller)) {};

public://Command      
	void solve(void);
	std::vector<double> calculate_error_norms(const Grid& grid) const;


private:
	void controll_time_step(const double current_time, double& time_step) const;

private:
	double current_time_ = 0.0;
	std::unique_ptr<Semi_Discrete_Equation> semi_discrete_equation_;
	std::unique_ptr<Time_Discrete_Scheme> time_discrete_scheme_;
	std::unique_ptr<Solve_End_Controller> end_controller_;
	std::unique_ptr<Solve_Post_Controller> post_controller_;
};