#include "../INC/Discrete_Equation.h"

void Discrete_Equation::solve(void)
{
    size_t current_iter = 0;
    double current_time = 0.0;
    Post_Processor::syncronize_solution_time(current_time);

    LOG << "================================================================================\n";
    LOG << "\t\t\t\t Solving\n";
    LOG << "================================================================================\n" << LOG.print_;

    SET_TIME_POINT;
    while (true)
    {
        if (this->end_controller_->is_time_to_end(current_iter, current_time))
        {            
            Post_Processor::post_solution();//post
            break;
        }

        if (this->post_controller_->is_time_to_post(current_iter, current_time))
        {
            Post_Processor::post_solution();//post
        }

        SET_TIME_POINT;
        auto time_step = this->semi_discrete_equation_->calculate_time_step();
        this->controll_time_step(current_time, time_step);

        this->time_discrete_scheme_->update(*this->semi_discrete_equation_, time_step);
        current_time += time_step;
        current_iter++;

        LOG << std::left << std::setw(5);
        LOG << "Iter:" << current_iter << "\t";
        LOG << "Time:" << current_time << "(" << this->end_controller_->progress_percentage_str(current_iter, current_time) << ")";
        LOG << "Computation cost:" << std::to_string(GET_TIME_DURATION) << "s \n" << LOG.print_;
    }

    LOG << "\n================================================================================\n";
    LOG << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
    LOG << "================================================================================\n\n" << LOG.print_;
    
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