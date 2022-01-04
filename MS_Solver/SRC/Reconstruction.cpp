#include "../INC/Reconstruction.h"

//#include "../INC/Post_Processor.h" //for debug
void Hierarchical_Limiting_DG::reconstruct(Discrete_Solution_DG& discrete_solution) const
{
	const auto num_cells = discrete_solution.num_cells();

	this->stability_criterion_->precaclulate(discrete_solution);
	this->indicator_->precalculate(discrete_solution);
	this->limiter_->precalculate(discrete_solution);

	////debug
	//std::vector<double> cell_type(num_cells);
	//for (uint cell_index = 0; cell_index < num_cells; ++cell_index)
	//{
	//	cell_type[cell_index] = static_cast<double>(this->indicator_->indicate(discrete_solution, cell_index, *this->stability_criterion_));
	//}

	//Post_Processor::record_solution();
	//Post_Processor::record_variables("cell_type", cell_type);
	//Post_Processor::post_solution();
	////

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