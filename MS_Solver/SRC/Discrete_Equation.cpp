#include "../INC/Discrete_Equation.h"

Discrete_Equation::Discrete_Equation(const Configuration& configuration)
{
    const auto time_discrete_scheme_name = configuration.get("time_discrete_scheme");
    this->time_discrete_scheme_ = Time_Discrete_Scheme_Factory::make(time_discrete_scheme_name);


}


void Discrete_Equation::solve(void)
{
    double current_time = 0.0;
    Post_Processing::syncronize_solution_time(current_time);

    Log::content_ << "================================================================================\n";
    Log::content_ << "\t\t\t\t Solving\n";
    Log::content_ << "================================================================================\n";

    SET_TIME_POINT;
    while (true) {
        if (this->solve_controller_->is_time_to_end(current_time)) {
            //Post_Processing::post_solution(solutions, "final");//post
            break;
        }

        //if (this->solve_controller_->is_time_to_post(current_time))
            //Post_Processing::post_solution(solutions);//post
            //Post_Processing::post_condition_ = true;

        SET_TIME_POINT;
        auto time_step = this->semi_discrete_equation_.calculate_time_step();
        this->solve_controller_->controll_time_step(current_time, time_step);

        this->time_discrete_scheme_->update(this->semi_discrete_equation_, time_step);
        current_time += time_step;

        Log::content_ << "computation cost: " << std::to_string(GET_TIME_DURATION) << "s \n";
        Log::print();
    }

    Log::content_ << "\n================================================================================\n";
    Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
    Log::content_ << "================================================================================\n\n";
    Log::print();
}