#include "../INC/Discrete_Equation.h"

Discrete_Equation::Discrete_Equation(const Configuration& configuration)
{
    this->end_controller_ = Solve_End_Controller_Factory::make(configuration);
    this->post_controller_ = Solve_Post_Controller_Factory::make(configuration);    
    this->time_discrete_scheme_ = Time_Discrete_Scheme_Factory::make(configuration);
}


void Discrete_Equation::solve(void)
{
    size_t current_iter = 0;
    double current_time = 0.0;
    Post_Processing::syncronize_solution_time(current_time);

    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Solving\n";
    Log::content_ << "================================================================================\n";

    SET_TIME_POINT;
    while (true) {
        if (this->end_controller_->is_time_to_end(current_iter, current_time)) {
            //Post_Processing::post_solution(solutions, "final");//post
            break;
        }

        //if (this->post_controller_->is_time_to_post(current_time))
            //Post_Processing::post_solution(solutions);//post
            //Post_Processing::post_condition_ = true;

        SET_TIME_POINT;
        auto time_step = this->semi_discrete_equation_.calculate_time_step();
        this->controll_time_step(current_time, time_step);

        this->time_discrete_scheme_->update(this->semi_discrete_equation_, time_step);
        current_time += time_step;
        current_iter++;

        Log::content_ << std::left << std::setw(5);
        Log::content_ << "Iter:" << current_iter << "\t";
        Log::content_ << "Time:" << current_time << "(" << this->end_controller_->progress_percentage_str(current_iter, current_time) << ")";
        Log::content_ << "Computation cost:" << std::to_string(GET_TIME_DURATION) << "s \n";
        Log::print();
    }

    Log::content_ << "\n================================================================================\n";
    Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
    Log::content_ << "================================================================================\n\n";
    Log::print();
}

void Discrete_Equation::controll_time_step(const double current_time, double& time_step) const
{
    if (this->end_controller_->is_need_to_controll_time_step(current_time, time_step)) 
    {
        this->end_controller_->controll_time_step(current_time, time_step);
        return;
    }

    if (this->post_controller_->is_need_to_controll_time_step(current_time, time_step))
        this->post_controller_->controll_time_step(current_time, time_step);
}