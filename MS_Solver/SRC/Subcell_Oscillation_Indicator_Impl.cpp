#include "../INC/Subcell_Oscillation_Indicator_Impl.h"

Subcell_Oscillation_Indicator_Default::Subcell_Oscillation_Indicator_Default(const Grid& grid,
	std::unique_ptr<Trouble_Boundary_Indicator>&& troubled_boundary_indicator,
	std::unique_ptr<Shock_Indicator>&& shock_indicator,
	std::unique_ptr<Discontinuity_Indicator>&& discontinuity_indicator)
	: num_cells_(grid.num_cells())
	, troubled_boundary_indicator_(std::move(troubled_boundary_indicator))
	, shock_indicator_(std::move(shock_indicator))
	, discontinuity_indicator_(std::move(discontinuity_indicator)) {};

void Subcell_Oscillation_Indicator_Default::check(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->cell_index_to_has_typeI_oscillation_table_.begin(), this->cell_index_to_has_typeI_oscillation_table_.end(), false);
	std::fill(this->cell_index_to_has_typeII_oscillation_table_.begin(), this->cell_index_to_has_typeII_oscillation_table_.end(), false);

	this->check_num_troubled_faces(discrete_solution);
	this->shock_indicator_->check(discrete_solution);
	this->discontinuity_indicator_->check(discrete_solution);
	
	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto num_troubled_faces = this->cell_index_to_num_troubled_faces_[cell_index];
		const auto near_shock = this->shock_indicator_->is_near_shock(cell_index);
		const auto near_discontinuity = this->discontinuity_indicator_->is_near_discontinuity(cell_index);

		if (this->typeI_threshold_number_ <= num_troubled_faces && near_shock)
		{
			this->cell_index_to_has_typeI_oscillation_table_[cell_index] = true;
		}

		if (this->typeII_threshold_number_ <= this->cell_index_to_num_troubled_faces_[cell_index] && near_discontinuity)
		{
			this->cell_index_to_has_typeII_oscillation_table_[cell_index] = true;
		}
	}
}

void Subcell_Oscillation_Indicator_Default::check_num_troubled_faces(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->cell_index_to_num_troubled_faces_.begin(), this->cell_index_to_num_troubled_faces_.end(), 0);

	const auto infc_index_to_scaled_avg_sol_jump_table = this->measurer_->measure_infc_index_to_scaled_average_solution_jump_table(discrete_solution);

	for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
	{
		const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto characteristic_length = this->infc_index_to_characteristic_length_table_[infc_index];
		const auto threshold_value = std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));

		const auto average_solution_jump = infc_index_to_scaled_avg_sol_jump_table[infc_index];
		if (threshold_value < average_solution_jump)
		{
			this->cell_index_to_num_troubled_faces_[oc_index]++;
			this->cell_index_to_num_troubled_faces_[nc_index]++;
		}
	}
}
