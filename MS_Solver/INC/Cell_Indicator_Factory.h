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
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_F_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_F_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {        
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_1_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_21_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_21_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_22_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_22_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_3_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_41_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_41_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_11_1_42_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_11_1_42_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_12_1_T_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_12_1_T_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_22_1_T_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        MLP_Indicator MLP_indicator(grid, discrete_solution, criterion_equation_index);
        auto subcell_oscillation_indicator = Subcell_Oscillation_Indicator_Factory::make_22_1_T_indicator(governing_equation_name, grid, discrete_solution, criterion_equation_index);

        return std::make_unique<hMLP_BD_Indicator>(std::move(MLP_indicator), std::move(subcell_oscillation_indicator));
    }
};