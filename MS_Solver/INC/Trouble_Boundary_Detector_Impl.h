#pragma once
#include "Indicator.h"
#include "Measuring_Function.h"

class Trouble_Boundary_Indicator_Type1 : public Trouble_Boundary_Indicator
{
public:
    Trouble_Boundary_Indicator_Type1(const Grid& grid, std::unique_ptr<Average_Solution_Jump_Measurer>&& measurer)
        : num_infcs_(grid.num_inner_faces())
        , infc_index_to_characteristic_length_table_(grid.inner_face_index_to_characteristic_length_table())
        , infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
        , measurer_(std::move(measurer))
    {
        this->cell_index_to_num_troubled_boundaries_table_.resize(grid.num_cells());
    }

public:
    void check(const Discrete_Solution_DG& discrete_solution)
    {
		std::fill(this->cell_index_to_num_troubled_boundaries_table_.begin(), this->cell_index_to_num_troubled_boundaries_table_.end(), 0);

		const auto infc_index_to_scaled_avg_sol_jump_table = this->measurer_->measure_infc_index_to_scaled_average_solution_jump_table(discrete_solution);

		for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
		{
			const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

			const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
			const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
			const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

			const auto characteristic_length = this->infc_index_to_characteristic_length_table_[infc_index];
			const auto threshold_value = std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));

			const auto scaled_avg_sol_jump = infc_index_to_scaled_avg_sol_jump_table[infc_index];
			if (threshold_value < scaled_avg_sol_jump)
			{
				this->cell_index_to_num_troubled_boundaries_table_[oc_index]++;
				this->cell_index_to_num_troubled_boundaries_table_[nc_index]++;
			}
		}
    }


private:
    uint num_infcs_ = 0;
    std::vector<double> infc_index_to_characteristic_length_table_;
    std::vector<std::pair<uint, uint>> infc_index_to_oc_nc_index_pair_table_;

    std::unique_ptr<Average_Solution_Jump_Measurer> measurer_;
};