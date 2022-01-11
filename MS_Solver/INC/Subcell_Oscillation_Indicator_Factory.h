#pragma once
#include "Discontinuity_Indicator_Factory.h"
#include "Shock_Indicator_Impl.h"
#include "Trouble_Boundary_Indicator_Factory.h"

//T_S_D
//T : trouble boundary indicator type code
//S : shock indicator type code
//D : discontinutiy indicator type code
class Subcell_Oscillation_Indicator_Factory
{
public:
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_T_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_true_indicator();

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_F_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_false_indicator();

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_1_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type1_indicator(grid,discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_21_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_21_indicator(grid,discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_22_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_22_indicator(grid, discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_3_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_type3_indicator(grid, discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_41_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_41_indicator(grid, discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_11_1_42_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_11_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_42_indicator(grid, discrete_solution);

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_12_1_T_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_12_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_true_indicator();

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
    static std::unique_ptr<Subcell_Oscillation_Indicator> make_22_1_T_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        auto trouble_boundary_indicator = Trouble_Boundary_Indicator_Factory::make_22_indicator(grid, discrete_solution, criterion_solution_index);
        auto shock_indicator = Shock_Indicator_Factory::make_type1_indicator(grid, discrete_solution, governing_equation_name);
        auto discontinuity_indicator = Discontinuity_Indicator_Factory::make_always_true_indicator();

        return std::make_unique<Subcell_Oscillation_Indicator>(grid, std::move(trouble_boundary_indicator), std::move(shock_indicator), std::move(discontinuity_indicator));
    }
};
