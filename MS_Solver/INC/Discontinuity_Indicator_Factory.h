#pragma once
#include "Discontinuity_Indicator_Impl.h"
#include "Measuring_Function_Impl.h"

class Discontinuity_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicator> make_always_true_indicator(void)
    {
        return std::make_unique<Always_True_Discontinuity_Indicator>();
    }
    static std::unique_ptr<Discontinuity_Indicator> make_always_false_indicator(void)
    {
        return std::make_unique<Always_False_Discontinuity_Indicator>();
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Type1_Discontinuity_Indicator>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type2_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Type2_Discontinuity_Indicator>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type3_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Type3_Discontinuity_Indicator>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type4_1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto density_index = 0;
        auto measurer = std::make_unique<Default_Average_Solution_Jump_Measurer>(grid, discrete_solution, density_index);
        return std::make_unique<Type4_Discontinuity_Indicator>(grid, discrete_solution, std::move(measurer));
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type4_2_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto density_index = 0;
        auto measurer = std::make_unique<Scaled_Average_Solution_Jump_Measurer>(grid, discrete_solution, density_index);
        return std::make_unique<Type4_Discontinuity_Indicator>(grid, discrete_solution, std::move(measurer));
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};