#pragma once

#include "Reconstruction.h"
#include "Limiter_Impl.h"
#include "Stability_Criterion_Impl.h"

#include "Indicator_Impl.h"
#include "Indicating_Function_Impl.h"

class Reconstruction_DG_Factory//static class
{
public://Query
    static std::unique_ptr<Reconstruction_DG> make_unique(const Configuration& configuration, const Grid& grid, Discrete_Solution_DG& discrete_solution);

private:
    Reconstruction_DG_Factory(void) = delete;
};