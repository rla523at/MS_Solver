#pragma once
#include "Indicator.h"

class Limiter
{
public://Command
    virtual void precalculate(const Discrete_Solution_DG& discrete_solution) abstract;

public://Query
    virtual void limit(const ushort cell_index, const Cell_Type Cell_Type, Discrete_Solution_DG& discrete_solution, ushort& projection_degree, const MLP_Criterion_Base& stability_criterion) const abstract;
    
    bool is_end(void) const { return this->is_end_; };

protected:
    mutable bool is_end_ = false;
};
