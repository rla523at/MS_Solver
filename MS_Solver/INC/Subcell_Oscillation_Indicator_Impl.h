#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Subcell_Oscillation_Indicator_Default : public Subcell_Oscillation_Indicator
{
public:
    Subcell_Oscillation_Indicator_Default(const Grid& grid,
        std::unique_ptr<Trouble_Boundary_Indicator>&& troubled_boundary_indicator,
        std::unique_ptr<Shock_Indicator>&& shock_indicator,
        std::unique_ptr<Discontinuity_Indicator>&& discontinuity_indicator);

public://Command
    void check(const Discrete_Solution_DG& discrete_solution);

private:
    static constexpr auto typeI_threshold_number_ = 2;
    static constexpr auto typeII_threshold_number_ = 1;
    
    uint num_cells_ = 0;
    
    std::unique_ptr<Trouble_Boundary_Indicator> troubled_boundary_indicator_;
    std::unique_ptr<Shock_Indicator> shock_indicator_;
    std::unique_ptr<Discontinuity_Indicator> discontinuity_indicator_;
};