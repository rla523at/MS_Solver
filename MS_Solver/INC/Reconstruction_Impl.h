#pragma once

#include "Reconstruction.h"
#include "Limiter_Impl.h"
#include "Stability_Criterion_Impl.h"

#include "Indicator_Impl.h"


#include "Indicating_Function_Impl.h"
#include "Measuring_Function_Impl.h"
class Discontinuity_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        if (ms::compare_icase(type_name, "Always_Flase"))
        {
            return std::make_unique<Always_False_Discontinuity_Indicator>();
        }
        else if (ms::compare_icase(type_name, "Always_True"))
        {
            return std::make_unique<Always_True_Discontinuity_Indicator>();
        }
        else if (ms::compare_icase(type_name, "Heuristic"))
        {
            const auto num_cells = grid.num_cells();
            constexpr auto allowable_diff = 0.1; // 10%까지 discontinuity가 아니라고 보겠다.
            
            std::vector<double> cell_index_to_threshold_value_table(num_cells, allowable_diff);
            auto measuring_function = std::make_unique<Max_Scaled_Difference_Measuring_Function>(grid, discrete_solution, criterion_solution_index);
            return std::make_unique<Discontinuity_Indicating_Function_Base>(cell_index_to_threshold_value_table, std::move(measuring_function));
        }
        else if (ms::compare_icase(type_name, "Extrapolation"))
        {
            //return std::make_unique<Extrapolation_Discontinuity_Indicator>(grid, discrete_solution, criterion_solution_index);
        }
        else
        {
            EXCEPTION(type_name + " is not supported discontinuity indicator type");
            return nullptr;
        }
    }

    static std::unique_ptr<Discontinuity_Indicating_Function> make_always_false(void)
    {
        return std::make_unique<Always_False_Discontinuity_Indicator>();
    }

    static std::unique_ptr<Discontinuity_Indicating_Function> make_always_true(void)
    {
        return std::make_unique<Always_True_Discontinuity_Indicator>();
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};

class Shock_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        if (ms::compare_icase(governing_equation_name, "Linear_Advection"))
        {
            return Discontinuity_Indicator_Factory::make_always_false();
        }
        else if (ms::compare_icase(governing_equation_name, "Burgers"))
        {
            constexpr auto criterion_solution_index = 0;
            return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, criterion_solution_index);
        }
        else if (ms::compare_icase(governing_equation_name, "Euler"))
        {
            const auto space_dimension = grid.space_dimension();

            if (space_dimension == 2)
            {
                return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_2D::pressure_index());
            }
            else if (space_dimension == 3)
            {
                return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, Euler_3D::pressure_index());
            }
            else
            {
                EXCEPTION("not supported space_dimension");
                return nullptr;
            }
        }
        else
        {
            EXCEPTION("not supported governing equation");
            return nullptr;
        }
    }

private:
    Shock_Indicator_Factory(void) = delete;
};

class Contact_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        constexpr auto solution_index = 0;
        return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, solution_index);
    }

private:
    Contact_Indicator_Factory(void) = delete;
};

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