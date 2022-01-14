#include "../INC/Trouble_Boundary_Indicator_Impl.h"

void Troubeld_Boudary_Indicator_with_Jump_Measurer_Base::check(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->cell_index_to_num_troubled_boundaries_table_.begin(), this->cell_index_to_num_troubled_boundaries_table_.end(), 0);

	const auto infc_index_to_sol_jump_table = this->face_jump_measurer_->measure_inner_face_index_to_solution_jump_table(discrete_solution);

	for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
	{
		const auto threshold_value = this->calculate_threshold_value(discrete_solution, infc_index);
		const auto sol_jump = infc_index_to_sol_jump_table[infc_index];

		if (threshold_value < sol_jump)
		{
			const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

			this->cell_index_to_num_troubled_boundaries_table_[oc_index]++;
			this->cell_index_to_num_troubled_boundaries_table_[nc_index]++;
		}
	}
}

double Trouble_Boundary_Indicator_Type1::calculate_threshold_value(const Discrete_Solution_DG& discrete_solution, const uint inner_face_index) const
{
	const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[inner_face_index];

	const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
	const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
	const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

	const auto characteristic_length = this->infc_index_to_characteristic_length_table_[inner_face_index];
	return std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));
}