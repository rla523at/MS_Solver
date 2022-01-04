#include "../INC/Indicator_Impl.h"

Cell_Type hMLP_BD_Indicator::indicate(const Discrete_Solution_DG& discrete_solution, const uint cell_index, const MLP_Criterion_Base& stability_criterion) const
{
    auto Cell_Type = this->MLP_indicator_.indicate(discrete_solution, cell_index, stability_criterion);

    switch (Cell_Type)
    {
    case Cell_Type::normal:
    {
        if (this->subcell_oscillation_indicator_.is_typeI_cell(cell_index) && this->shock_indicator_->has_discontinuity(cell_index))
        {
            Cell_Type = Cell_Type::typeI;
        }
        break;
    }
    case Cell_Type::smooth_extrema:
    {
        if (this->subcell_oscillation_indicator_.is_typeII_cell(cell_index) && this->contact_indicator_->has_discontinuity(cell_index))
        {
            Cell_Type = Cell_Type::typeII;
        }
        break;
    }
    default:
        break;
    }

    return Cell_Type;
}