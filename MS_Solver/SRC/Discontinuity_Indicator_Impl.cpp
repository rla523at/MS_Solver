#include "../INC/Discontinuity_Indicator_Impl.h"

Type1_Discontinuity_Indicator::Type1_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : measuring_function_(grid, discrete_solution, rho_index)
    , num_inner_faces_(grid.num_inner_faces())
    , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
{
    this->cell_index_to_near_discontinuity_table_.resize(grid.num_cells(), false);
};

void Type1_Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_near_discontinuity_table_.begin(), this->cell_index_to_near_discontinuity_table_.end(), false);

    const auto infc_index_to_scaled_average_difference_table = measuring_function_.measure_infc_index_to_scaled_average_difference_table(discrete_solution);
    for (uint infc_index = 0; infc_index < this->num_inner_faces_; ++infc_index)
    {
        const auto scaled_avg_diff = infc_index_to_scaled_average_difference_table[infc_index];

        if (this->threshold_number_ <= scaled_avg_diff)
        {
            const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];
            this->cell_index_to_near_discontinuity_table_[oc_index] = true;
            this->cell_index_to_near_discontinuity_table_[nc_index] = true;
        }
    }
};

Type2_Discontinuity_Indicator::Type2_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : num_cells_(grid.num_cells())
    , cell_index_to_characteristic_length_(grid.cell_index_to_characteristic_length_table())
    , measuring_function_(grid, discrete_solution, rho_index)
{
    this->cell_index_to_near_discontinuity_table_.resize(this->num_cells_, false);
}

void Type2_Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_near_discontinuity_table_.begin(), this->cell_index_to_near_discontinuity_table_.end(), false);

    const auto cell_index_to_extrapolation_differences = measuring_function_.measure_cell_index_to_extrapolation_differences(discrete_solution);
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto threshold_number = this->cell_index_to_characteristic_length_[cell_index];

        const auto& extrapolation_differences = cell_index_to_extrapolation_differences[cell_index];
        const auto max_extrapolation_diff = *std::max_element(extrapolation_differences.begin(), extrapolation_differences.end());

        if (threshold_number <= max_extrapolation_diff)
        {
            this->cell_index_to_near_discontinuity_table_[cell_index] = true;
        }
    }
};

Type3_Discontinuity_Indicator::Type3_Discontinuity_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : num_cells_(grid.num_cells())
    , measuring_function_(grid, discrete_solution)
{
    this->cell_index_to_near_discontinuity_table_.resize(this->num_cells_, false);
}

void Type3_Discontinuity_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_near_discontinuity_table_.begin(), this->cell_index_to_near_discontinuity_table_.end(), false);

    const auto cell_index_to_divergence_velocities_table = measuring_function_.measure_cell_index_to_divergence_velocities_table(discrete_solution);
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto& divergence_velocities = cell_index_to_divergence_velocities_table[cell_index];
        const auto max_divergence_velocity = *std::max_element(divergence_velocities.begin(), divergence_velocities.end());

        if (max_divergence_velocity <= this->threshold_number_)
        {
            this->cell_index_to_near_discontinuity_table_[cell_index] = true;
        }
    }
};

static std::unique_ptr<Discontinuity_Indicator> Discontinuity_Indicator_Factory::make_unique(const std::string& type_name, const Grid& grid, Discrete_Solution_DG& discrete_solution)
{
    if (ms::compare_icase(type_name, "type1"))
    {
        return std::make_unique<Type1_Discontinuity_Indicator>(grid, discrete_solution);
    }
    else if (ms::compare_icase(type_name, "type2"))
    {
        return std::make_unique<Type2_Discontinuity_Indicator>(grid, discrete_solution);
    }
    else if (ms::compare_icase(type_name, "type3"))
    {
        return std::make_unique<Type3_Discontinuity_Indicator>(grid, discrete_solution);
    }
    else
    {
        EXCEPTION("Not supported discontinuity indicator type");
        return nullptr;
    }
}