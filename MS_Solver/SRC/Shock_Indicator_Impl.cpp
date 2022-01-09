#include "../INC/Shock_Indicator_Impl.h"

Shock_Indicator_Type1::Shock_Indicator_Type1(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort pressure_index)
    : measuring_function_(grid, discrete_solution, pressure_index)
    , num_inner_faces_(grid.num_inner_faces())
    , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
{
    this->cell_index_to_near_shock_table_.resize(grid.num_cells(), false);
};

void Shock_Indicator_Type1::check(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_near_shock_table_.begin(), this->cell_index_to_near_shock_table_.end(), false);

    const auto infc_index_to_scaled_average_difference_table = measuring_function_.measure_infc_index_to_scaled_average_difference_table(discrete_solution);
    for (uint infc_index = 0; infc_index < this->num_inner_faces_; ++infc_index)
    {
        const auto scaled_avg_diff = infc_index_to_scaled_average_difference_table[infc_index];

        if (this->threshold_number_ <= scaled_avg_diff)
        {
            const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];
            this->cell_index_to_near_shock_table_[oc_index] = true;
            this->cell_index_to_near_shock_table_[nc_index] = true;
        }
    }
};

std::unique_ptr<Shock_Indicator> Shock_Indicator_Factory::make_unique(const std::string& governing_equation_name, const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
    const auto pressure_index = find_pressure_index(governing_equation_name, grid.space_dimension());

    if (pressure_index < 0)
    {
        return make_always_false_indicator();
    }

    if (ms::compare_icase(type_name, "type1"))
    {
        return std::make_unique<Shock_Indicator_Type1>(grid, discrete_solution, pressure_index);
    }
    else
    {
        EXCEPTION("not supported type");
        return nullptr;
    }
};