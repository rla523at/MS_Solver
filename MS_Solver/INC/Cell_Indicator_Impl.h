#pragma once
#include "Discontinuity_Indicator_Impl.h"
#include "Indicator.h"
#include "Indicating_Function.h"
#include "Shock_Indicator_Impl.h"

class MLP_Indicator : public Cell_Indicator
{
public:
    MLP_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_index);

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override;

public://Query
    Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override;

private:
    Cell_Type check_cell_type(const uint cell_index, const double* value_at_vertices, const MLP_Criterion_Base& criterion) const;
    bool is_constant(const double value, const double P0_value, const uint cell_index) const;
    bool is_satisfy_MLP_condition(const double P1_projected_value, const double allowable_min, const double allowable_max) const;
    bool is_smooth_extrema(const double value, const double higher_mode_value, const double P1_mode_value, const double allowable_min, const double allowable_max) const;

private:
    ushort criterion_equation_index_;
    uint num_cells_;

    std::vector<double> cell_index_to_volume_table_;
    std::vector<ushort> cell_index_to_num_vertices_table_;
    std::vector<std::vector<double>> cell_index_to_P1_projected_value_at_vertices_table_;

    //construction optimization
    static constexpr ushort num_max_vertices = 8;
    mutable std::array<double, num_max_vertices> value_at_vertices_;
    mutable std::array<double, num_max_vertices> P1_projected_value_at_vertices_;
};

class hMLP_BD_Indicator : public Cell_Indicator
{
public:
    hMLP_BD_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index
        , std::unique_ptr<Shock_Indicator>&& shock_indicator, std::unique_ptr<Discontinuity_Indicator>&& contact_indicator)
        : MLP_indicator_(grid, discrete_solution, criterion_equation_index)
        , subcell_oscillation_indicator_(grid, discrete_solution, criterion_equation_index)
        , shock_indicator_(std::move(shock_indicator))
        , discontinuity_indicator_(std::move(contact_indicator)) {};

public://Command
    void precalculate(const Discrete_Solution_DG& discrete_solution) override
    {
        this->MLP_indicator_.precalculate(discrete_solution);
        this->subcell_oscillation_indicator_.precalculate(discrete_solution);
        this->shock_indicator_->precalculate(discrete_solution);
        this->discontinuity_indicator_->precalculate(discrete_solution);
    }

public://Query
    Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const override;

private:
    MLP_Indicator MLP_indicator_;
    Subcell_Oscillation_Indicator subcell_oscillation_indicator_;
    std::unique_ptr<Shock_Indicator> shock_indicator_;
    std::unique_ptr<Discontinuity_Indicator> discontinuity_indicator_;
};

class Cell_Indicator_Factory//static class
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