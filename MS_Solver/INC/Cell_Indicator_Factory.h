#pragma once
#include "Cell_Indicator_Impl.h"
#include "Subcell_Oscillation_Indicator_Factory.h"

class Cell_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Cell_Indicator> make_hMLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        return std::make_unique<MLP_Indicator>(grid, discrete_solution, criterion_equation_index);
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_true_indicator();
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_Off_TypeII_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_false_indicator();
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type1_indicator(grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type2_indicator(grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type3_indicator(grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD4_1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_1_indicator(grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    //static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD4_2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    //{
    //    auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
    //    auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_2_indicator(grid, discrete_solution);
    //    return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    //}
};