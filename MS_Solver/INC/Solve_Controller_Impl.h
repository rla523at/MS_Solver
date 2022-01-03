#pragma once
#include "Configuration.h"
#include "Initial_Condition.h"
#include "Solve_Controller.h"

class Solve_End_Controller_By_Time : public Solve_End_Controller
{
public:
    Solve_End_Controller_By_Time(const double end_time) : end_time_(end_time) {};

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override;
    void controll_time_step(const double current_time, double& time_step) const override;
    bool is_time_to_end(const size_t current_iter, const double current_time) const override;
    std::string progress_percentage_str(const size_t current_iter, const double current_time) const override;

private:
    double end_time_ = 0.0;
};

class Solve_End_Controller_By_Iter : public Solve_End_Controller
{
public:
    Solve_End_Controller_By_Iter(const size_t end_iter) : end_iter_(end_iter) {};

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override;
    void controll_time_step(const double current_time, double& time_step) const override;
    bool is_time_to_end(const size_t current_iter, const double current_time) const override;
    std::string progress_percentage_str(const size_t current_iter, const double current_time) const override;

private:
    size_t end_iter_ = 0;
};

class Solve_End_Controller_Factory//static class
{
public:
    static std::unique_ptr<Solve_End_Controller> make_unique(const Configuration& configuration);

private:
    Solve_End_Controller_Factory(void) = delete;
};


class Solve_Post_Controller_Not_Use : public Solve_Post_Controller
{
public://Command
    void increase_num_post(void) override {};

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override { return false; };
    void controll_time_step(const double current_time, double& time_step) const override {};
    bool is_time_to_post(const size_t current_iter, const double current_time) const override { return false; };
    bool is_post_initial_solution(void) const override { return false; };
    bool is_post_final_solution(void) const override { return false; };
};

class Solve_Post_Controller_Inital_Final : public Solve_Post_Controller
{
public://Command
    void increase_num_post(void) override {};

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override { return false; };
    void controll_time_step(const double current_time, double& time_step) const override {};
    bool is_time_to_post(const size_t current_iter, const double current_time) const override { return false; };
};


class Solve_Post_Controller_By_Time : public Solve_Post_Controller
{
public:
    Solve_Post_Controller_By_Time(const double post_time_step);

public://Command
    void increase_num_post(void) override;

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override;
    void controll_time_step(const double current_time, double& time_step) const override;
    bool is_time_to_post(const size_t current_iter, const double current_time) const override;

private:
    void update_post_time(void);

private:
    double post_time_step_ = 0.0;
    double post_time_ = 0.0;
};

class Solve_Post_Controller_By_Iter : public Solve_Post_Controller
{
public:
    Solve_Post_Controller_By_Iter(const size_t post_iter_unit);

public://Command
    void increase_num_post(void) override;

public://Query
    bool is_need_to_controll_time_step(const double current_time, const double time_step) const override;
    void controll_time_step(const double current_time, double& time_step) const override;
    bool is_time_to_post(const size_t current_iter, const double current_time) const override;

private:
    void update_post_iter(void);

private:
    size_t post_iter_unit_ = 0;
    size_t post_iter_ = 0;
};

class Solve_Post_Controller_Factory//static class
{
public:
    static std::unique_ptr<Solve_Post_Controller> make_unique(const Configuration& configuration);

private:
    Solve_Post_Controller_Factory(void) = delete;
};