#pragma once
#include "Log.h"
#include "Exception.h"

#include <algorithm>
#include <type_traits>


class Solve_Controller abstract
{
public:
	virtual std::string end_condition_name(void) const abstract;
	virtual void controll_time_step(const double current_time, double& time_step) const abstract;
	virtual bool is_time_to_end(const double current_time) const abstract;
	virtual bool is_time_to_post(const double current_time) const abstract;
};


//
//
//enum class Controll_Condition
//{
//    by_time,
//    by_iter
//};
//
//
//class Solve_Controller
//{
//private:
//    using This_ = Solve_Controller;
//
//private:
//    inline static Controll_Condition end_condition_;
//    inline static Controll_Condition post_condition_;
//    inline static double end_condition_constant_;
//    inline static double post_condition_constant_;
//    inline static size_t num_post_ = 0;
//    inline static size_t num_iter_ = 0;
//
//private:
//    Solve_Controller(void) = delete;
//
//public:
//    static std::string end_condition_name(void) {
//        std::string str = "End at " + ms::double_to_string(end_condition_constant_);
//        if (end_condition_ == Controll_Condition::by_time)
//            return str + " time";
//        else
//            return str + " iter";
//    }
//
//    static void initialize(const Controll_Condition end_condition, const double end_condition_constant, const Controll_Condition post_condition, const double post_condition_constant) {
//        This_::end_condition_ = end_condition;
//        This_::end_condition_constant_ = end_condition_constant;
//        This_::post_condition_ = post_condition;
//        This_::post_condition_constant_ = post_condition_constant;
//        This_::num_post_ = 0;
//        This_::num_iter_ = 0;
//    };
//
//    static void controll_time_step(const double current_time, double& time_step) {
//        Log::content_ << "Iter:" << std::left << std::setw(5) << This_::num_iter_ << "\t";
//        Log::content_ << "current time: " << std::to_string(current_time) << " ";        
//
//        if (This_::end_condition_ == Controll_Condition::by_time)
//            Log::content_ << std::setprecision(2) << std::fixed << "(" << current_time * 100 / This_::end_condition_constant_ << "%)  \t" << std::defaultfloat << std::setprecision(6);
//        else
//            Log::content_ << std::setprecision(2) << std::fixed << "(" << This_::num_iter_ * 100 / This_::end_condition_constant_ << "%)  \t" << std::defaultfloat << std::setprecision(6);
//
//        Log::print();
//
//        This_::num_iter_++;
//
//        if (This_::end_condition_ == Controll_Condition::by_time) {
//            const double expect_time = current_time + time_step;
//            if (This_::end_condition_constant_ <= expect_time) {
//                const auto exceed_time = expect_time - This_::end_condition_constant_;
//                time_step -= exceed_time;
//
//                return;
//            }
//        }
//
//        if (This_::post_condition_ == Controll_Condition::by_time) {            
//            const auto post_time_step = This_::post_condition_constant_;
//
//            const auto target_time = (This_::num_post_ + 1) * post_time_step;
//            const auto expect_time = current_time + time_step;
//            if (target_time <= expect_time) {
//                const auto exceed_time = expect_time - target_time;
//                time_step -= exceed_time;
//                return;
//            }
//        }
//
//        //This_::num_iter_++;
//    }
//
//    static bool is_time_to_end(const double current_time) {                
//        if (This_::end_condition_ == Controll_Condition::by_time) {
//            if (ms::compare_double(current_time, This_::end_condition_constant_)) {
//                Log::content_ << "Iter:" << std::left << std::setw(5) << This_::num_iter_ << "\t";
//                Log::content_ << "current time: " << std::to_string(current_time) << " (100.00%)  \tend";
//                return true;
//            }
//            else
//                return false;
//        }
//        else {
//            if (This_::num_iter_ == This_::end_condition_constant_) {
//                Log::content_ << "Iter:" << std::left << std::setw(5) << This_::num_iter_ << "\t";
//                Log::content_ << "current time: " << std::to_string(current_time) << " (100.00%)  \tend";
//                return true;
//            }
//            else
//                return false;
//        }
//    }
//
//    static bool is_time_to_post(const double current_time) {
//        if (This_::post_condition_ == Controll_Condition::by_time) {
//            const auto post_time_step = This_::post_condition_constant_;
//            const auto target_time = (This_::num_post_ + 1) * post_time_step;
//
//            if (ms::compare_double(current_time, target_time)) {
//                This_::num_post_++;
//                return true;
//            }
//            else
//                return false;
//        }
//        else {
//            const auto target_iter = (This_::num_post_ + 1) * This_::post_condition_constant_;
//            if (This_::num_iter_ == target_iter) {
//                This_::num_post_++;
//                return true;
//            }
//            else
//                return false;
//        }
//    }
//};