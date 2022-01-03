#pragma once
#include <string>

class Solve_End_Controller
{
public://Query
    virtual bool is_need_to_controll_time_step(const double current_time, const double time_step) const abstract;
    virtual void controll_time_step(const double current_time, double& time_step) const abstract;
    virtual bool is_time_to_end(const size_t current_iter, const double current_time) const abstract;
    virtual std::string progress_percentage_str(const size_t current_iter, const double current_time) const abstract;
};

class Solve_Post_Controller
{
public://Command
    virtual void increase_num_post(void) abstract;

public://Query
    virtual bool is_need_to_controll_time_step(const double current_time, const double time_step) const abstract;
    virtual void controll_time_step(const double current_time, double& time_step) const abstract;
    virtual bool is_time_to_post(const size_t current_iter, const double current_time) const abstract;
    virtual bool is_post_initial_solution(void) const { return true; };
    virtual bool is_post_final_solution(void) const { return true; };

protected:
    size_t num_post_ = 0;
};