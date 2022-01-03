#pragma once

#include <cmath>
#include <vector>

#include "Stability_Criterion.h"

enum class cell_type
{
    normal = 0,
    smooth_extrema = 1,
    trouble = 2,
    typeI = 3,
    typeII = 4
};

class Indicator
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual cell_type indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const abstract;
};