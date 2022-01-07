#include "../INC/Indicator.h"

Subcell_Oscillation_Indicator::Subcell_Oscillation_Indicator(const Grid& grid, Discrete_Solution_DG& discrete_solution, const ushort criterion_equation_index)
	: num_infcs_(grid.num_inner_faces())
	, infc_index_to_characteristic_length_table_(grid.inner_face_index_to_characteristic_length_table())
	, infc_index_to_oc_nc_index_pair_table_(grid.inner_face_index_to_oc_nc_index_pair_table())
	, cell_index_to_num_troubled_faces_(grid.num_cells())
	, measuring_function_(grid, discrete_solution, criterion_equation_index) {};

void Subcell_Oscillation_Indicator::precalculate(const Discrete_Solution_DG& discrete_solution)
{
	std::fill(this->cell_index_to_num_troubled_faces_.begin(), this->cell_index_to_num_troubled_faces_.end(), 0);

	const auto infc_index_to_average_solution_jump_table = this->measuring_function_.measure_infc_index_to_average_solution_jump_table(discrete_solution);

	for (uint infc_index = 0; infc_index < this->num_infcs_; ++infc_index)
	{
		const auto [oc_index, nc_index] = this->infc_index_to_oc_nc_index_pair_table_[infc_index];

		const auto oc_solution_degree = discrete_solution.solution_degree(oc_index);
		const auto nc_solution_degree = discrete_solution.solution_degree(nc_index);
		const auto max_solution_degree = (std::max)(oc_solution_degree, nc_solution_degree);

		const auto characteristic_length = this->infc_index_to_characteristic_length_table_[infc_index];
		//const auto threshold_value = std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));

		//test
		const auto oc_P0 = discrete_solution.calculate_P0_nth_solution(oc_index, 0);
		const auto nc_P0 = discrete_solution.calculate_P0_nth_solution(nc_index, 0);
		const auto avg = 0.5 * (oc_P0 + nc_P0);
		const auto threshold_value = avg * std::pow(characteristic_length, 0.5 * (max_solution_degree + 1));
		//

		const auto average_solution_jump = infc_index_to_average_solution_jump_table[infc_index];
		if (threshold_value < average_solution_jump)
		{
			this->cell_index_to_num_troubled_faces_[oc_index]++;
			this->cell_index_to_num_troubled_faces_[nc_index]++;
		}
	}
}