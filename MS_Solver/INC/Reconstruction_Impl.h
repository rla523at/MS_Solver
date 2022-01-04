#pragma once

#include "Reconstruction.h"
#include "Limiter_Impl.h"
#include "Stability_Criterion_Impl.h"

#include "Indicator_Impl.h"
#include "Indicating_Function_Impl.h"

class Indicator_Factory
{
public:
    static std::unique_ptr<Indicator> make_hMLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        return std::make_unique<hMLP_Indicator>(grid, discrete_solution, criterion_equation_index);
    }
    static std::unique_ptr<Indicator> make_hMLP_BD_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Heuristic", grid, discrete_solution);
        auto contact_indicator = Discontinuity_Indicator_Factory::make_always_true();
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
    }
    static std::unique_ptr<Indicator> make_Improved_hMLP_BD1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Heuristic", grid, discrete_solution);
        auto contact_indicator = Contact_Indicator_Factory::make_unique("Heuristic", governing_equation_name, grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
    }
    static std::unique_ptr<Indicator> make_Improved_hMLP_BD2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Heuristic", grid, discrete_solution);
        auto contact_indicator = Contact_Indicator_Factory::make_unique("Extrapolation", governing_equation_name, grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
    }
};

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};