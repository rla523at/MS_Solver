#pragma once
#include "Indicating_Function.h"

class Not_Use_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Not_Use_Discontinuity_Indicator(const Grid& grid)
        : Discontinuity_Indicator(grid, NULL) {};

public:
    void precalculate(const Discrete_Solution_DG& discrete_solution) override {};
};

class Heuristic_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Heuristic_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
        : Discontinuity_Indicator(grid, criterion_solution_index) 
    {
        //for precalculation
        discrete_solution.precalculate_cell_P0_basis_values();
    };

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;
};

class Extrapolation_Discontinuity_Indicator : public Discontinuity_Indicator
{
public:
    Extrapolation_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

private:
    std::vector<double> cell_index_to_volume_reciprocal_table_;
    std::vector<double> cell_index_to_threshold_value_table_;
    std::vector<Euclidean_Vector> cell_index_to_QW_v_table_;
};

class Discontinuity_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index)
    {
        if (ms::compare_icase(type_name, "Not_Use"))
        {
            return std::make_unique<Not_Use_Discontinuity_Indicator>(grid);
        }
        else if (ms::compare_icase(type_name, "Heuristic"))
        {
            return std::make_unique<Heuristic_Discontinuity_Indicator>(grid, discrete_solution, criterion_solution_index);
        }
        else if (ms::compare_icase(type_name, "Extrapolation"))
        {
            return std::make_unique<Extrapolation_Discontinuity_Indicator>(grid, discrete_solution, criterion_solution_index);
        }
        else
        {
            EXCEPTION("not supported discontinuity indicator type");
        }
    }

private:
    Discontinuity_Indicator_Factory(void) = delete;
};

class Shock_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        if (ms::compare_icase(governing_equation_name, "Linear_Advection"))
        {
            return Discontinuity_Indicator_Factory::make_unique("Not_Use", grid, discrete_solution, NULL);
        }

        auto solution_index = 0;
        if (ms::compare_icase(governing_equation_name, "Euler"))
        {
            const auto space_dimension = grid.space_dimension();

            if (space_dimension == 2)
            {
                solution_index = Euler_2D::pressure_index();
            }
            else if (space_dimension == 3)
            {
                solution_index = Euler_3D::pressure_index();
            }
            else
            {
                EXCEPTION("not supported space_dimension");
            }
        }

        return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, solution_index);
    }

private:
    Shock_Indicator_Factory(void) = delete;
};

class Contact_Indicator_Factory
{
public:
    static std::unique_ptr<Discontinuity_Indicator> make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
    {
        if (ms::compare_icase(governing_equation_name, "Linear_Advection") || ms::compare_icase(governing_equation_name, "Burgers"))
        {
            return Discontinuity_Indicator_Factory::make_unique("Not_Use", grid, discrete_solution, NULL);
        }

        auto solution_index = 0;
        if (ms::compare_icase(governing_equation_name, "Euler"))
        {
            const auto space_dimension = grid.space_dimension();

            if (space_dimension == 2)
            {
                solution_index = Euler_2D::pressure_index();
            }
            else if (space_dimension == 3)
            {
                solution_index = Euler_3D::pressure_index();
            }
            else
            {
                EXCEPTION("not supported space_dimension");
            }
        }

        return Discontinuity_Indicator_Factory::make_unique(type_name, grid, discrete_solution, solution_index);
    }

private:
    Shock_Indicator_Factory(void) = delete;
};