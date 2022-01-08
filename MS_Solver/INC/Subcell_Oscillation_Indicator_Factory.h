#pragma once
#include "Discontinuity_Indicator_Factory.h"
#include "Shock_Indicator_Impl.h"
#include "Subcell_Oscillation_Indicator_Impl.h"

class Subcell_Oscillation_Indicator_Factory
{
public:
    std::unique_ptr<Subcell_Oscillation_Indicator> make_improved_type_4_1(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto measurer = std::make_unique<Default_Average_Solution_Jump_Measurer>(grid, discrete_solution, criterion_equation_index);
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_1_indicator(grid, discrete_solution);
        return std::make_unique<Subcell_Oscillation_Indicator_Default>(grid, std::move(measurer), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
	std::unique_ptr<Subcell_Oscillation_Indicator> make_improved_type_4_2(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto measurer = std::make_unique<Default_Average_Solution_Jump_Measurer>(grid, discrete_solution, criterion_equation_index);
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_2_indicator(grid, discrete_solution);
        return std::make_unique<Subcell_Oscillation_Indicator_Default>(grid, std::move(measurer), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
};
