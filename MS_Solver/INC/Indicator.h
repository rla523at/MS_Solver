#pragma once
#include <cmath>


class P1_Projected_MLP_Trouble_Vertex_Indicator
{
private:
    P1_Projected_MLP_Trouble_Vertex_Indicator(void) = delete;

public:
    static bool is_satisfy(const double P1_projected_value, const double allowable_min, const double allowable_max)
    {
        return allowable_min <= P1_projected_value && P1_projected_value <= allowable_max;
    }
};


class MLP_Smooth_Extrema_Vertex_Indicator
{
private:
    MLP_Smooth_Extrema_Vertex_Indicator(void) = delete;

public:
    static bool is_smooth_extrema(const double solution, const double higher_mode_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max)
    {
        if (P1_mode_solution > 0 && higher_mode_solution < 0 && solution > allowable_min)
        {
            return true;
        }
        else if (P1_mode_solution < 0 && higher_mode_solution > 0 && solution < allowable_max)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
};


class Constant_Region_Vertex_Indicator
{
private:
    Constant_Region_Vertex_Indicator(void) = delete;

public:
    static bool is_constant(const double solution, const double P0_solution, const double volume)
    {
        const auto constant_criterion = (std::max)(1.0E-3 * std::abs(P0_solution), volume);
        return std::abs(solution - P0_solution) <= constant_criterion;
    }
};