#pragma once
#include "Indicating_Function_Impl.h"

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
        auto contact_indicator = Contact_Indicator_Factory::make_unique("Hueristc", governing_equation_name, grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
    }
    static std::unique_ptr<Indicator> make_Improved_hMLP_BD2_Indicator(const std::string& governing_equation_name, const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
    {
        auto shock_indicator = Shock_Indicator_Factory::make_unique(governing_equation_name, "Heuristic", grid, discrete_solution);
        auto contact_indicator = Contact_Indicator_Factory::make_unique("Extrapolation", governing_equation_name, grid, discrete_solution);
        return std::make_unique<hMLP_BD_Indicator>(grid, discrete_solution, criterion_equation_index, std::move(shock_indicator), std::move(contact_indicator));
    }
};