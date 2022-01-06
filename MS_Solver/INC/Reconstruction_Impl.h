#pragma once

#include "Reconstruction.h"
#include "Limiter_Impl.h"
#include "Stability_Criterion_Impl.h"
#include "Indicating_Function_Impl.h"
#include "Measuring_Function_Impl.h"
#include "Cell_Indicator_Impl.h"

class Discontinuity_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_solution_index);
    static std::unique_ptr<Discontinuity_Indicating_Function> make_always_false(void);
    static std::unique_ptr<Discontinuity_Indicating_Function> make_always_true(void);

private:
    Discontinuity_Indicator_Factory(void) = delete;
};

class Shock_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Shock_Indicator_Factory(void) = delete;
};

class Contact_Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Discontinuity_Indicating_Function> make_unique(const std::string& type_name, const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Contact_Indicator_Factory(void) = delete;
};

class Indicator_Factory//static class
{
public://Query
    static std::unique_ptr<Cell_Indicator> make_hMLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_hMLP_BD_Off_TypeII_indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD1_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD3_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
    static std::unique_ptr<Cell_Indicator> make_Improved_hMLP_BD4_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index);
};

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};