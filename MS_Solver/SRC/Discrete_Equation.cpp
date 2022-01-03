#include "../INC/Discrete_Equation.h"

void Discrete_Equation::solve(void)
{
    size_t current_iter = 0;
    Post_Processor::syncronize_solution_time(this->current_time_);
    Post_Processor::post_grid();

    if (this->post_controller_->is_post_initial_solution())
    {
        Post_Processor::record_solution();
        Post_Processor::post_solution("initial");
    }

    while (true)
    {
        if (this->end_controller_->is_time_to_end(current_iter, this->current_time_))
        {
            if (this->post_controller_->is_post_final_solution())
            {
                Post_Processor::record_solution();
                Post_Processor::post_solution("final");
            }
            break;
        }

        if (this->post_controller_->is_time_to_post(current_iter, this->current_time_))
        {
            Post_Processor::record_solution();
            Post_Processor::post_solution();
            this->post_controller_->increase_num_post(); 
        }

        Profiler::set_time_point();

        auto time_step = this->semi_discrete_equation_->calculate_time_step();
        this->controll_time_step(this->current_time_, time_step);

        this->time_discrete_scheme_->update(*this->semi_discrete_equation_, time_step);
        this->current_time_ += time_step;
        current_iter++;

        LOG << std::left;
        LOG << std::setw(10) << "Iter:" + std::to_string(current_iter);
        LOG << std::setw(25) << "Time:" + std::to_string(this->current_time_) + "(" + this->end_controller_->progress_percentage_str(current_iter, this->current_time_) + ")";
        LOG << std::setw(25) << "Computation cost:" + std::to_string(Profiler::get_time_duration()) + "s \n" << LOG.print_;
    }    
}

void Discrete_Equation::controll_time_step(const double current_time, double& time_step) const
{
    if (this->end_controller_->is_need_to_controll_time_step(current_time, time_step)) 
    {
        this->end_controller_->controll_time_step(current_time, time_step);
        return;
    }

    if (this->post_controller_->is_need_to_controll_time_step(current_time, time_step))
    {
        this->post_controller_->controll_time_step(current_time, time_step);
    }
}

std::vector<double> Discrete_Equation::calculate_error_norms(const Grid& grid) const
{
    return this->semi_discrete_equation_->calculate_error_norms(grid, this->current_time_);
}