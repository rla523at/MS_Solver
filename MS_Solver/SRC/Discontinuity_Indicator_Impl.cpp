#include "../INC/Discontinuity_Indicator_Impl.h"

Discontinuity_Indicator_Type1::Discontinuity_Indicator_Type1(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : measuring_function_(grid, discrete_solution, rho_index)
    , num_inner_faces_(grid.num_inner_faces())
    , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
{
    this->cell_index_to_near_discontinuity_table_.resize(grid.num_cells(), false);
};

void Discontinuity_Indicator_Type1::check(const Discrete_Solution_DG& discrete_solution)
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

Discontinuity_Indicator_Type2::Discontinuity_Indicator_Type2(const Grid& grid, Discrete_Solution_DG& discrete_solution, std::unique_ptr<Extrapolation_Jump_Measurer>&& extrapolation_jump_measurer)
    : num_cells_(grid.num_cells())
    , cell_index_to_characteristic_length_(grid.cell_index_to_characteristic_length_table())
    , extrapolation_jump_measurer_(std::move(extrapolation_jump_measurer))
{
    this->cell_index_to_near_discontinuity_table_.resize(this->num_cells_, false);
}


//#include "../INC/Post_Processor.h"
void Discontinuity_Indicator_Type2::check(const Discrete_Solution_DG& discrete_solution)
{
    ////debug
    //std::vector<double> avgs(this->num_cells_);
    ////

    std::fill(this->cell_index_to_near_discontinuity_table_.begin(), this->cell_index_to_near_discontinuity_table_.end(), false);

    const auto cell_index_to_extrapolation_jumps = extrapolation_jump_measurer_->measure_cell_index_to_extrapolation_jumps(discrete_solution);
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto threshold_number = this->cell_index_to_characteristic_length_[cell_index];

        const auto& extrapolation_jumps = cell_index_to_extrapolation_jumps[cell_index];
        const auto num_face_share_cell = extrapolation_jumps.size();

        double sum = 0.0;
        for (ushort j = 0; j < num_face_share_cell; ++j)
        {
            sum += extrapolation_jumps[j];
        }
        const auto avg = sum / static_cast<double>(num_face_share_cell);

        ////debug
        //avgs[cell_index] = avg;
        ////

        if (threshold_number <= avg)
        {
            this->cell_index_to_near_discontinuity_table_[cell_index] = true;
        }
    }

    ////debug
    //Post_Processor::record_solution();
    //Post_Processor::record_variables("avgs", avgs);
    //Post_Processor::post_solution();
    ////
};

Discontinuity_Indicator_Type3::Discontinuity_Indicator_Type3(const Grid& grid, Discrete_Solution_DG& discrete_solution)
    : num_cells_(grid.num_cells())
    , measuring_function_(grid, discrete_solution)
{
    this->cell_index_to_near_discontinuity_table_.resize(this->num_cells_, false);
}

void Discontinuity_Indicator_Type3::check(const Discrete_Solution_DG& discrete_solution)
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

Discontinuity_Indicator_Type4::Discontinuity_Indicator_Type4(const Grid& grid, Discrete_Solution_DG& discrete_solution, std::unique_ptr<Face_Jump_Measurer>&& measurer)
    : num_cells_(grid.num_cells())
    , num_infcs_(grid.num_inner_faces())
    , cell_index_to_characteristic_length_table_(grid.cell_index_to_characteristic_length_table())
    , cell_index_to_num_infc_table_(grid.cell_index_to_num_inner_faces_table())
    , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
    , face_jump_measurer_(std::move(measurer))
{
    this->cell_index_to_near_discontinuity_table_.resize(this->num_cells_, false);
};

void Discontinuity_Indicator_Type4::check(const Discrete_Solution_DG& discrete_solution)
{
    std::fill(this->cell_index_to_near_discontinuity_table_.begin(), this->cell_index_to_near_discontinuity_table_.end(), false);

    std::vector<double> cell_index_to_sum_of_scaled_avg_sol_jump_table(this->num_cells_);

    const auto infc_index_to_scaled_avg_sol_jump_table = this->face_jump_measurer_->measure_inner_face_index_to_solution_jump_table(discrete_solution);
    for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
    {
        const auto scaled_avg_sol_jump = infc_index_to_scaled_avg_sol_jump_table[infc_index];

        const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];
        cell_index_to_sum_of_scaled_avg_sol_jump_table[oc_index] += scaled_avg_sol_jump;
        cell_index_to_sum_of_scaled_avg_sol_jump_table[nc_index] += scaled_avg_sol_jump;
    }    
    
    for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
    {
        const auto num_infcs = this->cell_index_to_num_infc_table_[cell_index];
        const auto sum = cell_index_to_sum_of_scaled_avg_sol_jump_table[cell_index];
        const auto average = sum / static_cast<double>(num_infcs);

        const auto threshold_number = this->cell_index_to_characteristic_length_table_[cell_index];

        if (threshold_number <= average)
        {
            this->cell_index_to_near_discontinuity_table_[cell_index] = true;
        }
    }
}