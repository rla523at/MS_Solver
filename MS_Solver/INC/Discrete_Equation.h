#pragma once
#include "Time_Discrete_Scheme.h"
#include "Post.h"
#include "Solve_Controller.h"

class Discrete_Equation
{
public:
    Discrete_Equation(const Configuration& configuration);

public://Command
    void solve(void);

private:
    void controll_time_step(const double current_time, double& time_step) const;

private:
    std::unique_ptr<Solve_End_Controller> end_controller_;
    std::unique_ptr<Solve_Post_Controller> post_controller_;

    std::unique_ptr<Time_Discrete_Scheme> time_discrete_scheme_;
    Semi_Discrete_Equation semi_discrete_equation_;
};


//#include "Debugger.h"
//#include "Semi_Discrete_Equation.h"
//#include "Time_Integral_Method.h"
//#include "Time_Step_Method.h"
//#include "Solve_Controller.h"
//#include "Post_Processing.h"

//template <typename Time_Integral_Method>
//class Discrete_Equation
//{
//private:
//    static_require(ms::is_time_integral_method<Time_Integral_Method>, "It should be time integral method");
//
//public: 
//    template<typename Time_Step_Method, typename Semi_Discrete_Equation, typename Solution>
//    static void solve(Semi_Discrete_Equation& semi_discrete_equation, std::vector<Solution>& solutions) {
        //double current_time = 0.0;
        //Post_Processing::syncronize_solution_time(current_time);

        //Log::content_ << "================================================================================\n";
        //Log::content_ << "\t\t\t\t Solving\n";
        //Log::content_ << "================================================================================\n";
        //                
        //Post_Processing::post_condition_ = true;
        //semi_discrete_equation.reconstruct(solutions);

        //SET_TIME_POINT;
        //while (true) {
        //    if (this->solve_controller_->is_time_to_end(current_time)) {
        //        Post_Processing::post_solution(solutions, "final");//post
        //        break;
        //    }                      

        //    
        //    if (this->solve_controller_->is_time_to_post(current_time)) 
        //        //Post_Processing::post_solution(solutions);//post
        //        Post_Processing::post_condition_ = true;

        //    SET_TIME_POINT;
        //    auto time_step = semi_discrete_equation.calculate_time_step<Time_Step_Method>(solutions);
        //    this->solve_controller_->controll_time_step(current_time, time_step);

        //    Time_Integral_Method::update_solutions(semi_discrete_equation, solutions, time_step);
        //    current_time += time_step;

        //    Log::content_ << "computation cost: " << std::to_string(GET_TIME_DURATION) << "s \n";
        //    Log::print();
        //}

        //Log::content_ << "\n================================================================================\n";
        //Log::content_ << "\t\t\t Total ellapsed time: " << GET_TIME_DURATION << "s\n";
        //Log::content_ << "================================================================================\n\n";
        //Log::print();
//    }
//};