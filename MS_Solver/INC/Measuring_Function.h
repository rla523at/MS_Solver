#pragma once
#include "Discrete_Solution.h"

class Discontinuity_Measuring_Function
{
public:
    virtual std::vector<double> measure(const Discrete_Solution_DG& discrete_solution) const abstract;
};
