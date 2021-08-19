#pragma once
#include "Log.h"

#include <algorithm>
#include <type_traits>


//class SEC {}; // Solve End Condition
//
//template<double target_time>
//class End_By_Time : public SEC
//{
//public:
//    static std::string name(void) { return "End at "+ std::to_string(target_time) + " Time"; };
//
//    static bool inspect(const double current_time, double& time_step) {
//        static count current_iter = 0;
//
        //Log::content_ << "current time: " << std::to_string(current_time) + "s  ";
        //Log::content_ << std::fixed << std::setprecision(3) << "(" << current_time * 100 / target_time << "%)\n" << std::defaultfloat << std::setprecision(6);
        //Log::content_ << "Iter:" << std::left << std::setw(5) << ++current_iter << "\t";
//
//        const double expect_time = current_time + time_step;
//        if (target_time <= expect_time) {
//            const auto exceed_time = expect_time - target_time;
//            time_step -= exceed_time;
//            current_iter = 0;
//            return true;
//        }
//        else
//            return false;
//    }
//
//    static double inspect_exceed_time_step(const double current_time, double& time_step) {        
//        const double expect_time = current_time + time_step;
//        if (target_time <= expect_time)
//            return expect_time - target_time;
//        else
//            return std::numeric_limits<double>::max();
//    }
//
//};
//
//template<count target_iter>
//class End_By_Iter : public SEC
//{
//public:
//    static std::string name(void) { return "End at " + std::to_string(target_iter) + " Iter"; };
//
//    static bool inspect(const double current_time, double& time_step) {
//        static count current_iter = 0;
//
//        Log::content_ << "current time: " << std::to_string(current_time) + "s  ";
//        Log::content_ << std::fixed << std::setprecision(3) << "(" << current_iter++ * 100 / target_iter << "%)\n" << std::defaultfloat << std::setprecision(6);
//        Log::content_ << "Iter:" << std::left << std::setw(5) << current_iter << "\t";
//                
//        if (current_iter == target_iter) {
//            current_iter = 0;
//            return true;
//        }
//        else
//            return false;
//    }
//};

//
//class SPC {};   // Solve Post Condition
//
//template<double post_time_step>
//class Post_By_Time : public SPC
//{
//public:
//    static bool inspect(const double current_time, double& time_step, ushort& num_post) {
//        const auto target_time = (num_post + 1) * post_time_step;
//        const auto expect_time = current_time + time_step;
//        if (target_time <= expect_time) {
//            const auto exceed_time = expect_time - target_time;
//            time_step -= exceed_time;
//            num_post++;
//            return true;
//        }
//        else
//            return false;
//    }
//};
//
//template<ushort post_iter>
//class Post_By_Iter : public SPC
//{
//public:
//    static bool inspect(const double current_time, double& time_step, ushort& num_post) {
//        static uint iter = 1;
//
//        if (iter++ == post_iter) {
//            iter = 1;
//            num_post++;
//            return true;
//        }
//        else
//            return false;
//    }
//};
//
//namespace ms {
//    template <typename T>
//    inline constexpr bool is_solve_end_condtion = std::is_base_of_v<SEC, T>;
//    template <typename T>
//    inline constexpr bool is_solve_post_condtion = std::is_base_of_v<SPC, T>;
//
//}



enum class controll_condition
{
    by_time,
    by_iter
};


class Solve_Controller
{
private:
    using This_ = Solve_Controller;

private:
    inline static controll_condition end_condition_;
    inline static controll_condition post_condition_;
    inline static double end_condition_constant_;
    inline static double post_condition_constant_;
    inline static size_t num_post_ = 0;
    inline static size_t num_iter_ = 0;

private:
    Solve_Controller(void) = delete;

public:
    static std::string solve_end_condition_name(void) {
        std::string str = "End at " + ms::double_to_string(end_condition_constant_);
        if (end_condition_ == controll_condition::by_time)
            return str + " time";
        else
            return str + " iter";
    }

    static void initialize(const controll_condition end_condition, const double end_condition_constant, const controll_condition post_condition, const double post_condition_constant) {
        This_::end_condition_ = end_condition;
        This_::end_condition_constant_ = end_condition_constant;
        This_::post_condition_ = post_condition;
        This_::post_condition_constant_ = post_condition_constant;
    };

    static void controll_time_step(const double current_time, double& time_step) {
        Log::content_ << "Iter:" << std::left << std::setw(5) << This_::num_iter_++ << "\t";
        Log::content_ << "current time: " << std::to_string(current_time) << " ";

        if (This_::end_condition_ == controll_condition::by_time)
            Log::content_ << std::setprecision(2) << std::fixed << "(" << current_time * 100 / This_::end_condition_constant_ << "%)  \t" << std::defaultfloat << std::setprecision(6);
        else
            Log::content_ << std::setprecision(2) << std::fixed << "(" << This_::num_iter_ * 100 / This_::end_condition_constant_ << "%)  \t" << std::defaultfloat << std::setprecision(6);

        std::vector<double> exceed_times;
        
        if (This_::end_condition_ == controll_condition::by_time) {
            const double expect_time = current_time + time_step;
            if (This_::end_condition_constant_ <= expect_time)
                exceed_times.push_back(expect_time - end_condition_constant_);            
        }

        if (This_::post_condition_ == controll_condition::by_time) {            

            const auto post_time_step = This_::post_condition_constant_;

            const auto target_time = (This_::num_post_ + 1) * post_time_step;
            const auto expect_time = current_time + time_step;
            if (target_time <= expect_time) 
                exceed_times.push_back(expect_time - target_time);
        }

        if (!exceed_times.empty()) {
            const auto min_exceed_time = *std::min_element(exceed_times.begin(), exceed_times.end());
            time_step -= min_exceed_time;
        }
    }

    static bool is_time_to_end(const double current_time) {                
        if (This_::end_condition_ == controll_condition::by_time) 
            return current_time == This_::end_condition_constant_;
        else 
            return This_::num_iter_ == This_::end_condition_constant_;
    }

    static bool is_time_to_post(const double current_time) {
        if (This_::post_condition_ == controll_condition::by_time) {
            const auto post_time_step = This_::post_condition_constant_;
            const auto target_time = (This_::num_post_ + 1) * post_time_step;

            if (current_time == target_time) {
                This_::num_post_++;
                return true;
            }
            else
                return false;
        }
        else {
            const auto target_iter = (This_::num_post_ + 1) * This_::post_condition_constant_;
            if (This_::num_iter_ == target_iter) {
                This_::num_post_++;
                return true;
            }
            else
                return false;
        }
    }
};