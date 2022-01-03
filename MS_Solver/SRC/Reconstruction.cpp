#include "../INC/Reconstruction.h"

void Hierarchical_Limiting_DG::reconstruct(Discrete_Solution_DG& discrete_solution) const
{
	const auto num_cells = discrete_solution.num_cells();

	this->stability_criterion_->precaclulate(discrete_solution);
	this->indicator_->precalculate(discrete_solution);
	this->limiter_->precalculate(discrete_solution);

	for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	{
		auto projection_degree = discrete_solution.solution_degree(cell_index);

		while (true)
		{
			const auto Cell_Type = this->indicator_->indicate(discrete_solution, cell_index, *this->stability_criterion_);
			this->limiter_->limit(cell_index, Cell_Type, discrete_solution, projection_degree, *this->stability_criterion_);

			if (this->limiter_->is_end())
			{
				break;
			}
		}
	}
}

//void Test_Reconstuction_DG::reconstruct(Discrete_Solution_DG& discrete_solution)
//{
//	this->discontinuity_indicator_.precalculate(discrete_solution);
//	const auto& discontinuity_factor = this->discontinuity_indicator_.get_discontinuity_factor();
//	Post_Processor::record_solution();
//	Post_Processor::record_variables("discontinuity_factor", discontinuity_factor);
//	Post_Processor::post_solution();
//};