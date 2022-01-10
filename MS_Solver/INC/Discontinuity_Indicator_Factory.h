#pragma once
#include "Discontinuity_Indicator_Impl.h"
#include "Measuring_Function_Impl.h"

//D(F)
//D : Discontinuity Indicator Type
//(F) : Face Jump Measurer Type 
class Discontinuity_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicator> make_always_true_indicator(void)
    {
        return std::make_unique<Discontinuity_Indicator_Always_True>();
    }
    static std::unique_ptr<Discontinuity_Indicator> make_always_false_indicator(void)
    {
        return std::make_unique<Discontinuity_Indicator_Always_False>();
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Discontinuity_Indicator_Type1>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type2_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Discontinuity_Indicator_Type2>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type3_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        return std::make_unique<Discontinuity_Indicator_Type3>(grid, discrete_solution);
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type4_1_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto density_index = 0;
        auto measurer = std::make_unique<Face_Jump_Measurer_Type1>(grid, discrete_solution, density_index);
        return std::make_unique<Discontinuity_Indicator_Type4>(grid, discrete_solution, std::move(measurer));
    }
    static std::unique_ptr<Discontinuity_Indicator> make_type4_2_indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto density_index = 0;
        auto measurer = std::make_unique<Face_Jump_Measurer_Type2>(grid, discrete_solution, density_index);
        return std::make_unique<Discontinuity_Indicator_Type4>(grid, discrete_solution, std::move(measurer));
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};