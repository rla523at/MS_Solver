#include "../INC/Indicator.h"

Subcell_Oscillation_Indicator::Subcell_Oscillation_Indicator(const Grid& grid,
	std::unique_ptr<Trouble_Boundary_Indicator>&& troubled_boundary_indicator,
	std::unique_ptr<Shock_Indicator>&& shock_indicator,
	std::unique_ptr<Discontinuity_Indicator>&& discontinuity_indicator)
	: num_cells_(grid.num_cells())
	, trouble_boundary_indicator_(std::move(troubled_boundary_indicator))
	, shock_indicator_(std::move(shock_indicator))
	, discontinuity_indicator_(std::move(discontinuity_indicator)) {};

void Subcell_Oscillation_Indicator::check(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->cell_index_to_has_typeI_oscillation_table_.begin(), this->cell_index_to_has_typeI_oscillation_table_.end(), false);
	std::fill(this->cell_index_to_has_typeII_oscillation_table_.begin(), this->cell_index_to_has_typeII_oscillation_table_.end(), false);

	this->trouble_boundary_indicator_->check(discrete_solution);
	this->shock_indicator_->check(discrete_solution);
	this->discontinuity_indicator_->check(discrete_solution);

	for (uint cell_index = 0; cell_index < this->num_cells_; ++cell_index)
	{
		const auto num_troubled_bdrys = this->trouble_boundary_indicator_->num_troubled_boundary(cell_index);
		const auto near_shock = this->shock_indicator_->is_near_shock(cell_index);
		const auto near_discontinuity = this->discontinuity_indicator_->is_near_discontinuity(cell_index);

		if (this->typeI_threshold_number_ <= num_troubled_bdrys && near_shock)
		{
			this->cell_index_to_has_typeI_oscillation_table_[cell_index] = true;
		}

		if (this->typeII_threshold_number_ <= num_troubled_bdrys && near_discontinuity)
		{
			this->cell_index_to_has_typeII_oscillation_table_[cell_index] = true;
		}
	}
}
