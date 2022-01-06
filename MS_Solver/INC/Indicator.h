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

class Cell_Indicator
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual Cell_Type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const abstract;
};

class Discontinuity_Indicator
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual bool near_discontinuity(const uint cell_index) const abstract;

protected:
    std::vector<bool> cell_index_to_near_discontinuity_table_;
};

class Shock_Indicator
{
public:
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public:
    virtual bool near_shock(const uint cell_index) const abstract;

protected:
    std::vector<bool> cell_index_to_near_discontinuity_table_;
};

