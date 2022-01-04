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
    Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override
    {
        return this->MLP_indicator_.indicate(discrete_solution, cell_index, stability_criterion);
    }

private:    
    MLP_Indicator MLP_indicator_;
};

class hMLP_BD_Indicator : public Indicator
{
public:
    hMLP_BD_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index
        , std::unique_ptr<Discontinuity_Indicator>&& shock_indicator, std::unique_ptr<Discontinuity_Indicator>&& contact_indicator)
        : MLP_indicator_(grid, discrete_solution, criterion_equation_index)
        , subcell_oscillation_indicator_(grid, discrete_solution, criterion_equation_index)
        , shock_indicator_(std::move(shock_indicator))
        , contact_indicator_(std::move(contact_indicator)) {};

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override
    {
        this->MLP_indicator_.precalculate(discrete_solution);
        this->subcell_oscillation_indicator_.precalculate(discrete_solution);
        this->shock_indicator_->precalculate(discrete_solution);
        this->contact_indicator_->precalculate(discrete_solution);
    }

public://Query
    Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override;

private:
    MLP_Indicator MLP_indicator_;
    Subcell_Oscillation_Indicator subcell_oscillation_indicator_;
    std::unique_ptr<Discontinuity_Indicator> shock_indicator_;
    std::unique_ptr<Discontinuity_Indicator> contact_indicator_;
};