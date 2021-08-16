#pragma once
#include "Log.h"

#include <type_traits>


using count =  unsigned long long;

class SEC {}; // Solve End Condition

template<double target_time>
class End_By_Time : public SEC
{
public:
    static std::string name(void) { return "End at "+ std::to_string(target_time) + " Time"; };

    static bool inspect(const double current_time, double& time_step) {
        static count current_iter = 0;

        Log::content_ << "current time: " << std::to_string(current_time) + "s  ";
        Log::content_ << std::fixed << std::setprecision(3) << "(" << current_time * 100 / target_time << "%)\n" << std::defaultfloat << std::setprecision(6);
        Log::content_ << "Iter:" << std::left << std::setw(5) << ++current_iter << "\t";

        const double expect_time = current_time + time_step;
        if (target_time <= expect_time) {
            const auto exceed_time = expect_time - target_time;
            time_step -= exceed_time;
            return true;
        }
        else
            return false;
    }
};

template<count target_iter>
class End_By_Iter : public SEC
{
public:
    static std::string name(void) { return "End at " + std::to_string(target_iter) + " Iter"; };

    static bool inspect(const double current_time, double& time_step) {
        static count current_iter = 0;

        Log::content_ << "current time: " << std::to_string(current_time) + "s  ";
        Log::content_ << std::fixed << std::setprecision(3) << "(" << current_iter++ * 100 / target_iter << "%)\n" << std::defaultfloat << std::setprecision(6);
        Log::content_ << "Iter:" << std::left << std::setw(5) << current_iter << "\t";
                
        if (current_iter == target_iter) 
            return true;
        else
            return false;
    }
};


class SPC {};   // Solve Post Condition

template<double post_time_step>
class Post_By_Time : public SPC
{
public:
    static bool inspect(const double current_time, double& time_step, ushort& num_post) {
        const auto target_time = (num_post + 1) * post_time_step;
        const auto expect_time = current_time + time_step;
        if (target_time <= expect_time) {
            const auto exceed_time = expect_time - target_time;
            time_step -= exceed_time;
            num_post++;
            return true;
        }
        else
            return false;
    }
};

template<ushort post_iter>
class Post_By_Iter : public SPC
{
public:
    static bool inspect(const double current_time, double& time_step, ushort& num_post) {
        static uint iter = 1;

        if (iter++ == post_iter) {
            iter = 1;
            num_post++;
            return true;
        }
        else
            return false;
    }
};

namespace ms {
    template <typename T>
    inline constexpr bool is_solve_end_condtion = std::is_base_of_v<SEC, T>;
    template <typename T>
    inline constexpr bool is_solve_post_condtion = std::is_base_of_v<SPC, T>;

}
