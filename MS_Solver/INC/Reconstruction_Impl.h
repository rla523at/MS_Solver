#pragma once
#include "Cell_Indicator_Impl.h"
#include "Limiter_Impl.h"
#include "Reconstruction.h"
#include "Stability_Criterion_Impl.h"

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};