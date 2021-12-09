#include "../INC/Solve_Controller.h"

bool Solve_End_Controller_By_Time::is_need_to_controll_time_step(const double current_time, const double time_step) const
{
    const double expect_time = current_time + time_step;
    return this->end_time_ < expect_time;
}

bool Solve_End_Controller_By_Time::is_time_to_end(const size_t current_iter, const double current_time) const
{
    return current_time == this->end_time_;
}

void Solve_End_Controller_By_Time::controll_time_step(const double current_time, double& time_step) const
{
    time_step = this->end_time_ - current_time;
}

std::string Solve_End_Controller_By_Time::progress_percentage_str(const size_t current_iter, const double current_time) const
{
    std::ostringstream oss;
    oss << std::setprecision(2) << std::fixed << current_time * 100 / this->end_time_ << "%";
    return oss.str();
}


bool Solve_End_Controller_By_Iter::is_need_to_controll_time_step(const double current_time, const double time_step) const
{
    return false;
}

bool Solve_End_Controller_By_Iter::is_time_to_end(const size_t current_iter, const double current_time) const
{
    return current_iter == this->end_iter_;
}

void Solve_End_Controller_By_Iter::controll_time_step(const double current_time, double& time_step) const
{
    EXCEPTION("controller by iter never controll time step");
    return;
}

std::string Solve_End_Controller_By_Iter::progress_percentage_str(const size_t current_iter, const double current_time) const
{
    std::ostringstream oss;
    oss << std::setprecision(2) << std::fixed << current_iter * 100 / this->end_iter_ << "%";
    return oss.str();
}

Solve_Post_Controller_By_Time::Solve_Post_Controller_By_Time(const double post_time_step)
    : post_time_step_(post_time_step) {
    this->update_post_time();
};


void Solve_Post_Controller_By_Time::increase_num_post(void)
{
    this->num_post_++;
    this->update_post_time();
}

bool Solve_Post_Controller_By_Time::is_need_to_controll_time_step(const double current_time, const double time_step) const
{
    const auto expect_time = current_time + time_step;
    return this->post_time_ < expect_time;
}

bool Solve_Post_Controller_By_Time::is_time_to_post(const size_t current_iter, const double current_time) const
{
    return current_time == this->post_time_;
}

void Solve_Post_Controller_By_Time::controll_time_step(const double current_time, double& time_step) const
{
    time_step = this->post_time_ - current_time;
}

void Solve_Post_Controller_By_Time::update_post_time(void)
{
    this->post_time_ = (this->num_post_ + 1) * this->post_time_step_;
}

Solve_Post_Controller_By_Iter::Solve_Post_Controller_By_Iter(const size_t post_iter_unit) 
    : post_iter_unit_(post_iter_unit) 
{
    this->update_post_iter();
};

void Solve_Post_Controller_By_Iter::increase_num_post(void)
{
    this->num_post_++;
    this->update_post_iter();
}

bool Solve_Post_Controller_By_Iter::is_need_to_controll_time_step(const double current_time, const double time_step) const
{
    return false;
}

bool Solve_Post_Controller_By_Iter::is_time_to_post(const size_t current_iter, const double current_time) const
{
    return current_iter == this->post_iter_;
}

void Solve_Post_Controller_By_Iter::controll_time_step(const double current_time, double& time_step) const
{
    EXCEPTION("controller by iter never controll time step");
    return;
}

void Solve_Post_Controller_By_Iter::update_post_iter(void)
{
    this->post_iter_ = (this->num_post_ + 1) * this->post_iter_unit_;
}

std::unique_ptr<Solve_End_Controller> Solve_End_Controller_Factory::make_unique(const Configuration& configuration)
{
    const auto& type_name = configuration.get_solve_end_controller_type();
    if (ms::contains_icase(type_name, "Time"))
    {
        const auto solve_end_time = configuration.solve_end_time();
        return std::make_unique<Solve_End_Controller_By_Time>(solve_end_time);
    }
    else if (ms::contains_icase(type_name, "Iter"))
    {
        const auto solve_end_iter = configuration.solve_end_iter();
        return std::make_unique<Solve_End_Controller_By_Iter>(solve_end_iter);
    }
    else
    {
        EXCEPTION("solve_end_controller_type in configuration file does not supported");
        return nullptr;
    }
}

std::unique_ptr<Solve_Post_Controller> Solve_Post_Controller_Factory::make_unique(const Configuration& configuration)
{
    const auto& type_name = configuration.get_solve_post_controller_type();
    if (ms::contains_icase(type_name, "Time"))
    {
        const auto post_time_step = configuration.solve_post_time_step();
        return std::make_unique<Solve_Post_Controller_By_Time>(post_time_step);
    }
    else if (ms::contains_icase(type_name, "Iter"))
    {
        const auto post_iter_unit = configuration.solve_post_iter_unit();
        return std::make_unique<Solve_Post_Controller_By_Iter>(post_iter_unit);
    }
    else
    {
        EXCEPTION("solve_end_controller_type in configuration file does not supported");
        return nullptr;
    }
}