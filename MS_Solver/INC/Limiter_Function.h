#pragma once
#include <vector>

using ushort = unsigned short;

class MLP_u1_Limiter
{
private:
    MLP_u1_Limiter(void) = delete;

public:
    static double limiter_function(const double P0_value, const std::vector<double>& P1_projected_values, const std::vector<std::pair<double, double>>& set_of_allowable_min_max)
    {
        const auto num_vertices = P1_projected_values.size();

        double limiting_value = 1.0;
        for (ushort j = 0; j < num_vertices; ++j)
        {
            const auto [allowable_min, allowable_max] = set_of_allowable_min_max[j];

            const auto P1_mode_criterion_value = P1_projected_values[j] - P0_value;

            limiting_value = (std::min)(limiting_value, MLP_u1_Limiter::calculate_limiting_value(P0_value, P1_mode_criterion_value, allowable_min, allowable_max));
        }

        return limiting_value;
    };

    static double calculate_limiting_value(const double P0_solution, const double P1_mode_solution, const double allowable_min, const double allowable_max) {
        if (P1_mode_solution == 0)
        {
            return 1.0;
        }

        if (P1_mode_solution < 0)
        {
            return (std::min)((allowable_min - P0_solution) / P1_mode_solution, 1.0);
        }
        else
        {
            return (std::min)((allowable_max - P0_solution) / P1_mode_solution, 1.0);
        }
    };
};
