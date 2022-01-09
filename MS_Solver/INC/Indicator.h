#pragma once
#include "Stability_Criterion.h"

enum class Cell_Type
{
    normal = 0,
    smooth_extrema = 1,
    trouble = 2,
    typeI = 3,
    typeII = 4
};

class Cell_Indicator abstract
{
public://Command
    virtual void check(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const abstract;
};

class Discontinuity_Indicator abstract
{
public://Command
    virtual void check(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual bool is_near_discontinuity(const uint cell_index) const { return cell_index_to_near_discontinuity_table_[cell_index]; };

protected:
    std::vector<bool> cell_index_to_near_discontinuity_table_;
};

class Shock_Indicator
{
public:
    virtual void check(const Discrete_Solution_DG& discrete_solution) abstract;

public:
    virtual bool is_near_shock(const uint cell_index) const { return cell_index_to_near_shock_table_[cell_index]; };

protected:
    std::vector<bool> cell_index_to_near_shock_table_;
};

class Trouble_Boundary_Indicator abstract
{
public:
    virtual void check(const Discrete_Solution_DG& discrete_solution) abstract;

public:
    ushort num_troubled_boundary(const uint cell_index)  const { return cell_index_to_num_troubled_boundaries_table_[cell_index]; };

protected:
    std::vector<ushort> cell_index_to_num_troubled_boundaries_table_;
};

class Subcell_Oscillation_Indicator
{
public:
    Subcell_Oscillation_Indicator(const Grid& grid,
        std::unique_ptr<Trouble_Boundary_Indicator>&& trouble_boundary_indicator,
        std::unique_ptr<Shock_Indicator>&& shock_indicator,
        std::unique_ptr<Discontinuity_Indicator>&& discontinuity_indicator);

public://Command
    void check(const Discrete_Solution_DG& discrete_solution);

public://Query
    bool is_typeI_cell(const uint cell_index) const { cell_index_to_has_typeI_oscillation_table_[cell_index]; };
    bool is_typeII_cell(const uint cell_index) const { cell_index_to_has_typeII_oscillation_table_[cell_index]; };

private:
    static constexpr auto typeI_threshold_number_ = 2;
    static constexpr auto typeII_threshold_number_ = 1;

    uint num_cells_ = 0;
    std::vector<bool> cell_index_to_has_typeI_oscillation_table_;
    std::vector<bool> cell_index_to_has_typeII_oscillation_table_;

    std::unique_ptr<Trouble_Boundary_Indicator> trouble_boundary_indicator_;
    std::unique_ptr<Shock_Indicator> shock_indicator_;
    std::unique_ptr<Discontinuity_Indicator> discontinuity_indicator_;

};