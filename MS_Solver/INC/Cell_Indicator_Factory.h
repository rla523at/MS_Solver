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
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_T_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_F_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_F_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {        
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_1_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_2_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_3_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_41_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_41_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD_11_1_42_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_42_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }

    //static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    //{
    //    auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
    //    auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type3_indicator(grid, discrete_solution);
    //    return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    //}
    //static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD4_1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    //{
    //    auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
    //    auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_1_indicator(grid, discrete_solution);
    //    return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    //}
    //static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD4_2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    //{
    //    auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Type1", grid, discrete_solution);
    //    auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type4_2_indicator(grid, discrete_solution);
    //    return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(discontinuity_indicator));
    //}
};