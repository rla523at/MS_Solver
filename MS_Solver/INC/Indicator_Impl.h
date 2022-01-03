#pragma once
#include "Indicating_Function.h"

class hMLP_Indicator : public Indicator
{
public:
    hMLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
        :MLP_indicator_(grid, discrete_solution, criterion_equation_index) {};

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override
    {
        this->MLP_indicator_.precalculate(discrete_solution);
    }

public://Query
    cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override
    {
        return this->MLP_indicator_.indicate(discrete_solution, cell_index, stability_criterion);
    }

private:    
    MLP_Indicator MLP_indicator_;
};

class hMLP_BD_Indicator : public Indicator
{
public:
    hMLP_BD_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index, const std::string& governing_equation_name)
        : MLP_indicator_(grid, discrete_solution, criterion_equation_index)
        , subcell_oscillation_indicator_(grid, discrete_solution, criterion_equation_index)
        , shock_indicator_(grid, discrete_solution, governing_equation_name) {};

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override
    {
        this->MLP_indicator_.precalculate(discrete_solution);
        this->subcell_oscillation_indicator_.precalculate(discrete_solution);
        this->shock_indicator_.precalculate(discrete_solution);
    }

public://Query
    cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override
    {
        auto cell_type = this->MLP_indicator_.indicate(discrete_solution, cell_index, stability_criterion);

        switch (cell_type)
        {
        case cell_type::normal:
        {
            if (this->subcell_oscillation_indicator_.is_typeI_cell(cell_index) && shock_indicator_.is_shock(cell_index))
            {
                cell_type = cell_type::typeI;
            }
            break;
        }
        case cell_type::smooth_extrema:
        {
            if (this->subcell_oscillation_indicator_.is_typeII_cell(cell_index))
            {
                cell_type = cell_type::typeII;
            }
            break;
        }
        default:
            break;
        }

        return cell_type;
    }

private:
    MLP_Indicator MLP_indicator_;
    Subcell_Oscillation_Indicator subcell_oscillation_indicator_;
    Shock_Indicator shock_indicator_;
};